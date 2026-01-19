#include "eastr/streaming_spurious_detector.hpp"
#include "eastr/self_aligner.hpp"
#include "eastr/bowtie2_runner.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>

namespace eastr {

struct StreamingSpuriousDetector::Impl {
    FastaIndex& fasta;
    std::string fasta_path;
    std::string bt2_index;
    ScoringParams scoring;
    AlgorithmParams params;
    std::unordered_map<std::string, int64_t> chrom_sizes;

    Impl(FastaIndex& f, const std::string& bt2, const ScoringParams& s, const AlgorithmParams& p)
        : fasta(f), fasta_path(f.get_path()), bt2_index(bt2), scoring(s), params(p) {
        chrom_sizes = fasta.get_chrom_sizes();
    }

    // Check if alignment spans both anchors (two-anchor alignment)
    static bool is_two_anchor_alignment(const AlignmentHit& hit, int overhang, int anchor) {
        bool c1 = (hit.r_en < overhang + anchor - 1);
        bool c2 = (hit.r_st > overhang - anchor);
        bool c3 = (hit.q_en < overhang + anchor - 1);
        bool c4 = (hit.q_st > overhang - anchor);
        bool c5 = (std::abs(hit.r_st - hit.q_st) > anchor * 2);

        return !(c1 || c2 || c3 || c4 || c5);
    }

    // Calculate linear distance between two sequences
    static int linear_distance(const std::string& s1, const std::string& s2) {
        if (s1.size() != s2.size()) {
            return -1;
        }

        int distance = 0;
        for (size_t i = 0; i < s1.size(); ++i) {
            if (s1[i] != s2[i]) {
                distance++;
            }
        }
        return distance;
    }

    // Check if a junction is spurious based on alignment results
    bool is_spurious_alignment(const JunctionKey& key,
                               const JunctionData& data,
                               const SequenceCache& seqs) const {
        if (!data.hit) {
            return false;  // No self-alignment, not spurious
        }

        const AlignmentHit& hit = *data.hit;
        int overhang = params.overhang;
        int anchor = params.anchor;
        int bt2_k = params.bt2_k;
        int min_dup_len = params.min_duplicate_exon_length;

        bool is_two_anchor = is_two_anchor_alignment(hit, overhang, anchor);

        // Get sequences
        auto rseq_it = seqs.find(data.jstart_region);
        auto qseq_it = seqs.find(data.jend_region);

        if (rseq_it == seqs.end() || qseq_it == seqs.end()) {
            return true;  // Can't verify, mark as spurious
        }

        const std::string& rseq = rseq_it->second;
        const std::string& qseq = qseq_it->second;

        if (is_two_anchor) {
            // Two-anchor alignment logic
            if ((data.seq1_count == 1 || data.seq2_count == 1) && data.seqh_count == 0) {
                return false;  // Unique alignment
            }
        } else {
            // One-anchor alignment logic

            // Check for duplicated exon
            if (hit.q_st - hit.r_st >= min_dup_len) {
                if (data.seqh_count == 0) {
                    return false;  // Duplicated exon, not spurious
                }
            }

            // Check alignment uniqueness
            if ((data.seq1_count < bt2_k) || (data.seq2_count < bt2_k)) {
                if (data.seqh_count == 0) {
                    return false;  // Partial alignment - no hybrid sequence found
                }

                // Check anchor region similarity
                if (static_cast<int>(rseq.size()) >= overhang + anchor &&
                    static_cast<int>(qseq.size()) >= overhang + anchor) {

                    std::string e5 = rseq.substr(overhang - anchor, anchor);
                    std::string i5 = rseq.substr(overhang, anchor);
                    std::string i3 = qseq.substr(overhang - anchor, anchor);
                    std::string e3 = qseq.substr(overhang, anchor);

                    int dist_e5i3 = linear_distance(e5, i3);
                    int dist_i5e3 = linear_distance(i5, e3);

                    if (dist_e5i3 >= 0 && dist_i5e3 >= 0) {
                        if (dist_e5i3 + dist_i5e3 > 2) {
                            return false;  // Anchors not similar enough - NOT spurious
                        }
                    }
                }
            }
        }

        return true;  // Default: spurious
    }
};

StreamingSpuriousDetector::StreamingSpuriousDetector(
    FastaIndex& fasta,
    const std::string& bt2_index,
    const ScoringParams& scoring,
    const AlgorithmParams& params,
    int batch_size,
    int num_threads)
    : pImpl_(std::make_unique<Impl>(fasta, bt2_index, scoring, params)),
      batch_size_(batch_size),
      num_threads_(num_threads) {
}

StreamingSpuriousDetector::~StreamingSpuriousDetector() = default;

void StreamingSpuriousDetector::add_junction(JunctionKey key, JunctionData data) {
    batch_.emplace_back(std::move(key), std::move(data));
}

std::vector<std::pair<JunctionKey, JunctionData>>
StreamingSpuriousDetector::process_batch(bool is_bam_input) {
    if (batch_.empty()) {
        return {};
    }

    std::vector<std::pair<JunctionKey, JunctionData>> spurious;

    // Convert batch to JunctionMap for compatibility with existing functions
    JunctionMap batch_map;
    for (auto& [key, data] : batch_) {
        batch_map[key] = std::move(data);
    }

    // Get flanking sequences for this batch (parallel for multiple threads)
    SequenceCache seqs = get_flanking_subsequences_parallel(
        batch_map, pImpl_->chrom_sizes, pImpl_->params.overhang,
        pImpl_->fasta_path, num_threads_);

    // Convert back to vector for parallel processing
    std::vector<std::pair<JunctionKey, JunctionData>> batch_vec;
    batch_vec.reserve(batch_map.size());
    for (auto& [key, data] : batch_map) {
        batch_vec.emplace_back(key, std::move(data));
    }

    // Self-align introns (parallel)
    auto self_aligned = get_self_aligned_introns_batch(
        batch_vec, seqs, pImpl_->params.overhang,
        pImpl_->scoring, pImpl_->params, num_threads_);

    // Separate low-score junctions from ones needing bowtie2
    JunctionMap introns_to_align;
    std::vector<std::pair<JunctionKey, JunctionData>> low_score_spurious;

    for (auto& [key, data] : self_aligned) {
        if (is_bam_input && data.total_score <= pImpl_->params.min_junc_score) {
            // Low-score junctions are automatically spurious
            low_score_spurious.emplace_back(key, std::move(data));
        } else {
            introns_to_align[key] = std::move(data);
        }
    }

    // Add low-score spurious to results
    for (auto& item : low_score_spurious) {
        spurious.push_back(std::move(item));
    }

    // Run bowtie2 alignment on remaining introns
    if (!introns_to_align.empty()) {
        Bowtie2Runner bt2_runner(pImpl_->bt2_index, num_threads_, pImpl_->params.bt2_k + 1);
        bt2_runner.align_probes(introns_to_align, seqs, pImpl_->params.overhang);
    }

    // Classify each junction
    for (const auto& [key, data] : introns_to_align) {
        if (pImpl_->is_spurious_alignment(key, data, seqs)) {
            spurious.emplace_back(key, data);
        }
    }

    // Update statistics
    total_processed_ += batch_.size();
    total_spurious_ += spurious.size();

    // Clear the batch
    batch_.clear();

    return spurious;
}

std::vector<std::pair<JunctionKey, JunctionData>>
StreamingSpuriousDetector::flush(bool is_bam_input) {
    return process_batch(is_bam_input);
}

} // namespace eastr
