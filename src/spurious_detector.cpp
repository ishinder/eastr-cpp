#include "eastr/spurious_detector.hpp"
#include "eastr/bed_parser.hpp"

#include <algorithm>
#include <cmath>

namespace eastr {

SpuriousDetector::SpuriousDetector(const ScoringParams& scoring,
                                   const AlgorithmParams& params,
                                   const std::string& bt2_index,
                                   const std::string& ref_fasta,
                                   int num_threads)
    : scoring_(scoring),
      params_(params),
      bt2_index_(bt2_index),
      fasta_(std::make_unique<FastaIndex>(ref_fasta)),
      num_threads_(num_threads) {
}

SpuriousDetector::~SpuriousDetector() = default;

bool SpuriousDetector::is_two_anchor_alignment(const AlignmentHit& hit,
                                                int overhang,
                                                int anchor) {
    // From get_spurious_introns.py:102-111
    bool c1 = (hit.r_en < overhang + anchor - 1);
    bool c2 = (hit.r_st > overhang - anchor);
    bool c3 = (hit.q_en < overhang + anchor - 1);
    bool c4 = (hit.q_st > overhang - anchor);
    bool c5 = (std::abs(hit.r_st - hit.q_st) > anchor * 2);

    return !(c1 || c2 || c3 || c4 || c5);
}

int SpuriousDetector::linear_distance(const std::string& s1, const std::string& s2) {
    if (s1.size() != s2.size()) {
        return -1;  // Invalid comparison
    }

    int distance = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i]) {
            distance++;
        }
    }
    return distance;
}

bool SpuriousDetector::is_spurious_alignment(const JunctionKey& key,
                                              const JunctionData& data,
                                              const SequenceCache& seqs) const {
    if (!data.hit) {
        return false;  // No self-alignment, not spurious
    }

    const AlignmentHit& hit = *data.hit;
    int overhang = params_.overhang;
    int anchor = params_.anchor;
    int bt2_k = params_.bt2_k;
    int min_dup_len = params_.min_duplicate_exon_length;

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
        // Two-anchor alignment logic (get_spurious_introns.py:118-124)
        // If either seq1 or seq2 has unique alignment AND hybrid unmapped -> NOT spurious
        if ((data.seq1_count == 1 || data.seq2_count == 1) && data.seqh_count == 0) {
            return false;  // Unique alignment
        }
    } else {
        // One-anchor alignment logic (get_spurious_introns.py:127-157)

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

            // Check anchor region similarity (when seqh_count > 0)
            // From Python get_spurious_introns.py:142-153
            if (static_cast<int>(rseq.size()) >= overhang + anchor &&
                static_cast<int>(qseq.size()) >= overhang + anchor) {

                // Extract anchor sequences
                std::string e5 = rseq.substr(overhang - anchor, anchor);  // Exon 5' end
                std::string i5 = rseq.substr(overhang, anchor);            // Intron 5' start
                std::string i3 = qseq.substr(overhang - anchor, anchor);  // Intron 3' end
                std::string e3 = qseq.substr(overhang, anchor);            // Exon 3' start

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

std::unordered_set<JunctionKey, JunctionKeyHash>
SpuriousDetector::load_trusted_junctions(const std::string& bed_path) {
    auto junctions = BedParser::parse(bed_path);

    std::unordered_set<JunctionKey, JunctionKeyHash> trusted;
    for (const auto& [key, _] : junctions) {
        trusted.insert(key);
    }

    return trusted;
}

JunctionMap SpuriousDetector::detect_spurious(
    const JunctionMap& junctions,
    bool is_bam_input,
    const std::optional<std::string>& trusted_bed_path) {

    // Load trusted junctions if provided
    std::unordered_set<JunctionKey, JunctionKeyHash> trusted;
    if (trusted_bed_path) {
        trusted = load_trusted_junctions(*trusted_bed_path);
    }

    // Filter out trusted junctions
    JunctionMap filtered;
    for (const auto& [key, data] : junctions) {
        if (trusted.find(key) == trusted.end()) {
            filtered[key] = data;
        }
    }

    // Get chromosome sizes
    auto chrom_sizes = fasta_->get_chrom_sizes();

    // Get flanking sequences (parallel)
    SequenceCache seqs = get_flanking_subsequences_parallel(
        filtered, chrom_sizes, params_.overhang, fasta_->get_path(), num_threads_);

    // Find self-aligned introns (parallel)
    // Convert to vector for batch processing
    std::vector<std::pair<JunctionKey, JunctionData>> batch;
    batch.reserve(filtered.size());
    for (auto& [key, data] : filtered) {
        batch.emplace_back(key, std::move(data));
    }

    auto aligned_batch = get_self_aligned_introns_batch(
        batch, seqs, params_.overhang, scoring_, params_, num_threads_);

    // Convert back to map
    JunctionMap self_aligned;
    for (auto& [key, data] : aligned_batch) {
        self_aligned[key] = std::move(data);
    }

    // Separate introns that need bowtie2 validation from low-score ones
    JunctionMap introns_to_align;
    JunctionMap spurious;

    for (auto& [key, data] : self_aligned) {
        if (is_bam_input && data.total_score <= params_.min_junc_score) {
            // Low-score junctions are automatically spurious
            spurious[key] = data;
        } else {
            introns_to_align[key] = data;
        }
    }

    // Run bowtie2 alignment on remaining introns
    if (!introns_to_align.empty()) {
        Bowtie2Runner bt2_runner(bt2_index_, num_threads_, params_.bt2_k + 1);
        bt2_runner.align_probes(introns_to_align, seqs, params_.overhang);
    }

    // Classify each junction
    for (const auto& [key, data] : introns_to_align) {
        if (is_spurious_alignment(key, data, seqs)) {
            spurious[key] = data;
        }
    }

    // Sort by position (for consistent output)
    return spurious;
}

} // namespace eastr
