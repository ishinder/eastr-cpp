#include "eastr/self_aligner.hpp"
#include "eastr/thread_pool.hpp"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <mutex>
#include <vector>

#ifdef EASTR_USE_MINIMAP2
extern "C" {
#include "minimap.h"
}
#endif

namespace eastr {

struct SelfAligner::Impl {
    ScoringParams scoring;
    AlgorithmParams params;

#ifdef EASTR_USE_MINIMAP2
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
#endif

    Impl(const ScoringParams& s, const AlgorithmParams& p) : scoring(s), params(p) {
#ifdef EASTR_USE_MINIMAP2
        // Initialize minimap2 options (similar to mappy defaults)
        mm_set_opt(0, &iopt, &mopt);

        // Set k-mer and window size from params
        iopt.k = static_cast<short>(params.kmer_size);
        iopt.w = static_cast<short>(params.window_size);

        // Set scoring parameters
        // mappy scoring format: [match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2]
        mopt.a = scoring.match_score;
        mopt.b = scoring.mismatch_penalty;
        mopt.q = scoring.gap_open1;
        mopt.e = scoring.gap_extend1;
        mopt.q2 = scoring.gap_open2;
        mopt.e2 = scoring.gap_extend2;

        // Set minimum chain score
        mopt.min_chain_score = params.min_chain_score;

        // We want only the best alignment
        mopt.best_n = 1;

        // Enable CIGAR generation
        mopt.flag |= MM_F_CIGAR;
#endif
    }

#ifdef EASTR_USE_MINIMAP2
    // Get score for aligning two bases (fallback for non-minimap2 builds)
#else
    int score(char a, char b) const {
        if (a == 'N' || a == 'n' || b == 'N' || b == 'n') {
            return -scoring.ambiguous_score;
        }
        if (a == b || (std::toupper(a) == std::toupper(b))) {
            return scoring.match_score;
        }
        return -scoring.mismatch_penalty;
    }
#endif
};

SelfAligner::SelfAligner(const ScoringParams& scoring, const AlgorithmParams& params)
    : pImpl_(std::make_unique<Impl>(scoring, params)) {
}

SelfAligner::~SelfAligner() = default;

SelfAligner::SelfAligner(SelfAligner&&) noexcept = default;
SelfAligner& SelfAligner::operator=(SelfAligner&&) noexcept = default;

std::optional<AlignmentHit> SelfAligner::align(const std::string& ref_seq,
                                                const std::string& query_seq) {
    if (ref_seq.empty() || query_seq.empty()) {
        return std::nullopt;
    }

#ifdef EASTR_USE_MINIMAP2
    // Build index from reference sequence
    const char* seqs[1] = { ref_seq.c_str() };
    const char* names[1] = { "ref" };

    mm_idx_t* mi = mm_idx_str(pImpl_->iopt.w, pImpl_->iopt.k, 0,
                               pImpl_->iopt.bucket_bits, 1, seqs, names);
    if (!mi) {
        return std::nullopt;
    }

    // Update mapping options with index info
    mm_mapopt_update(&pImpl_->mopt, mi);

    // CRITICAL: Match mappy behavior - don't filter high-occurrence seeds
    // mappy.pyx line 159 sets this for sequence-based indexes
    // Without this, minimap2 uses a low mid_occ threshold (calculated from index size)
    // which filters out many seeds and leads to different (often worse) alignments
    pImpl_->mopt.mid_occ = 1000;

    // Create thread buffer
    mm_tbuf_t* tbuf = mm_tbuf_init();

    // Map query against reference
    int n_regs = 0;
    mm_reg1_t* regs = mm_map(mi, static_cast<int>(query_seq.size()),
                              query_seq.c_str(), &n_regs, tbuf,
                              &pImpl_->mopt, nullptr);

    std::optional<AlignmentHit> result;

    // Find the first primary, forward-strand alignment
    for (int i = 0; i < n_regs; ++i) {
        mm_reg1_t* r = &regs[i];

        // Skip reverse strand alignments (like Python does)
        if (r->rev) {
            continue;
        }

        // Skip non-primary alignments
        if (r->parent != r->id) {
            continue;
        }

        // Build result
        AlignmentHit hit;
        hit.strand = 1;  // Forward
        hit.is_primary = true;

        // Positions (0-based)
        hit.r_st = r->rs;
        hit.r_en = r->re;
        hit.q_st = r->qs;
        hit.q_en = r->qe;

        // Match length
        hit.mlen = r->mlen;

        // Build CIGAR string and calculate score from CIGAR
        // Score formula: matches * match_score - mismatches * mismatch_penalty - gap_penalty
        // Gap penalty: min(gap_open1 + len*gap_extend1, gap_open2 + len*gap_extend2)
        // This matches minimap2's DP formula (not Python's buggy (len-1) formula)
        if (r->p && r->p->n_cigar > 0) {
            std::string cigar_str;
            int matches = 0;
            int mismatches = 0;
            int gap_penalty = 0;
            int r_pos = r->rs;  // Position in reference
            int q_pos = r->qs;  // Position in query

            for (uint32_t j = 0; j < r->p->n_cigar; ++j) {
                uint32_t c = r->p->cigar[j];
                int len = c >> 4;
                int op = c & 0xf;
                cigar_str += std::to_string(len);
                cigar_str += "MIDNSHP=XB"[op];

                if (op == 0) {  // M (match/mismatch)
                    // Count actual matches vs mismatches by comparing sequences
                    for (int k = 0; k < len; ++k) {
                        char r_base = std::toupper(ref_seq[r_pos + k]);
                        char q_base = std::toupper(query_seq[q_pos + k]);
                        if (r_base == q_base) {
                            ++matches;
                        } else {
                            ++mismatches;
                        }
                    }
                    r_pos += len;
                    q_pos += len;
                } else if (op == 1) {  // I (insertion in query)
                    // Gap penalty: min of two affine gap costs
                    int pen1 = pImpl_->scoring.gap_open1 + len * pImpl_->scoring.gap_extend1;
                    int pen2 = pImpl_->scoring.gap_open2 + len * pImpl_->scoring.gap_extend2;
                    gap_penalty += std::min(pen1, pen2);
                    q_pos += len;
                } else if (op == 2) {  // D (deletion in query)
                    // Gap penalty: min of two affine gap costs
                    int pen1 = pImpl_->scoring.gap_open1 + len * pImpl_->scoring.gap_extend1;
                    int pen2 = pImpl_->scoring.gap_open2 + len * pImpl_->scoring.gap_extend2;
                    gap_penalty += std::min(pen1, pen2);
                    r_pos += len;
                }
                // Skip N (intron), S (soft clip), H (hard clip), P (padding)
            }
            hit.cigar = cigar_str;

            // Calculate score using the correct formula
            hit.score = matches * pImpl_->scoring.match_score
                      - mismatches * pImpl_->scoring.mismatch_penalty
                      - gap_penalty;
        } else {
            // Fallback if no CIGAR available
            hit.score = (r->p) ? r->p->dp_max : r->score;
        }

        result = hit;
        break;  // Take first valid hit
    }

    // Cleanup
    for (int i = 0; i < n_regs; ++i) {
        if (regs[i].p) free(regs[i].p);
    }
    free(regs);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);

    return result;

#else
    // Fallback: Smith-Waterman local alignment
    const int m = static_cast<int>(ref_seq.size());
    const int n = static_cast<int>(query_seq.size());

    std::vector<std::vector<int>> H(m + 1, std::vector<int>(n + 1, 0));
    int max_score = 0;
    int max_i = 0, max_j = 0;

    const int gap_open = pImpl_->scoring.gap_open1;
    const int gap_extend = pImpl_->scoring.gap_extend1;

    std::vector<std::vector<int>> E(m + 1, std::vector<int>(n + 1, -1000000));
    std::vector<std::vector<int>> F(m + 1, std::vector<int>(n + 1, -1000000));

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match = H[i-1][j-1] + pImpl_->score(ref_seq[i-1], query_seq[j-1]);
            E[i][j] = std::max(E[i][j-1] - gap_extend, H[i][j-1] - gap_open - gap_extend);
            F[i][j] = std::max(F[i-1][j] - gap_extend, H[i-1][j] - gap_open - gap_extend);
            H[i][j] = std::max({0, match, E[i][j], F[i][j]});

            if (H[i][j] > max_score) {
                max_score = H[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    if (max_score < pImpl_->params.min_chain_score) {
        return std::nullopt;
    }

    // Traceback
    int i = max_i, j = max_j;
    std::string cigar_ops;

    while (i > 0 && j > 0 && H[i][j] > 0) {
        int current = H[i][j];
        int diag = H[i-1][j-1] + pImpl_->score(ref_seq[i-1], query_seq[j-1]);

        if (current == diag) {
            cigar_ops = 'M' + cigar_ops;
            --i; --j;
        } else if (current == E[i][j]) {
            cigar_ops = 'I' + cigar_ops;
            --j;
        } else if (current == F[i][j]) {
            cigar_ops = 'D' + cigar_ops;
            --i;
        } else {
            break;
        }
    }

    int r_st = i, r_en = max_i;
    int q_st = j, q_en = max_j;

    // Length filter
    if (std::min(r_en - r_st, q_en - q_st) < 15) {
        return std::nullopt;
    }

    AlignmentHit hit;
    hit.score = max_score;
    hit.strand = 1;
    hit.is_primary = true;
    hit.r_st = r_st;
    hit.r_en = r_en;
    hit.q_st = q_st;
    hit.q_en = q_en;

    // Compress CIGAR
    if (!cigar_ops.empty()) {
        std::string compressed;
        char prev_op = cigar_ops[0];
        int count = 1;
        for (size_t k = 1; k < cigar_ops.size(); ++k) {
            if (cigar_ops[k] == prev_op) {
                ++count;
            } else {
                compressed += std::to_string(count) + prev_op;
                prev_op = cigar_ops[k];
                count = 1;
            }
        }
        compressed += std::to_string(count) + prev_op;
        hit.cigar = compressed;
    }

    hit.mlen = std::min(hit.r_en - hit.r_st, hit.q_en - hit.q_st);
    return hit;
#endif
}

int SelfAligner::calc_alignment_score(const AlignmentHit& hit) const {
    // Use minimap2's raw alignment score directly.
    // Note: minimap2's r->score is computed from the actual DP alignment
    // using the scoring parameters we set (match, mismatch, gap penalties).
    //
    // Python's mappy recalculates from cs tag because mappy's hit.mlen
    // is the actual number of matching bases from the alignment,
    // while mm_reg1_t's mlen is "seeded exact match length" (from minimizers).
    //
    // Using the raw score is correct and equivalent to Python's recalculation
    // since minimap2 already computed it with our scoring parameters.
    return hit.score;
}

JunctionMap get_self_aligned_introns(
    const JunctionMap& introns,
    const SequenceCache& seqs,
    int overhang,
    const ScoringParams& scoring,
    const AlgorithmParams& params) {

    SelfAligner aligner(scoring, params);
    JunctionMap self_aligned;

    for (const auto& [key, data] : introns) {
        auto rseq_it = seqs.find(data.jstart_region);
        auto qseq_it = seqs.find(data.jend_region);

        if (rseq_it == seqs.end() || qseq_it == seqs.end()) {
            continue;
        }

        const std::string& rseq = rseq_it->second;
        const std::string& qseq = qseq_it->second;

        auto hit = aligner.align(rseq, qseq);

        if (!hit) {
            continue;
        }

        // Filter based on intron length and overlap
        int64_t intron_len = key.length();
        // Match Python's elif structure exactly:
        // - Condition 1: hit.r_st == intron_len
        // - Condition 2 (elif): overhang > intron_len AND hit.r_st - hit.q_st == intron_len
        // - Condition 3 (elif): hit.r_st >= overhang AND hit.q_en <= overhang
        // Condition 3 is only checked when overhang <= intron_len
        if (overhang * 2 >= intron_len) {
            if (hit->r_st == intron_len) {
                continue;
            } else if (overhang > intron_len) {
                if (hit->r_st - hit->q_st == intron_len) {
                    continue;
                }
            } else if (hit->r_st >= overhang && hit->q_en <= overhang) {
                continue;
            }
        }

        JunctionData new_data = data;
        new_data.hit = hit;
        self_aligned[key] = new_data;
    }

    return self_aligned;
}

std::vector<std::pair<JunctionKey, JunctionData>> get_self_aligned_introns_batch(
    std::vector<std::pair<JunctionKey, JunctionData>>& batch,
    const SequenceCache& seqs,
    int overhang,
    const ScoringParams& scoring,
    const AlgorithmParams& params,
    int num_threads) {

    std::vector<std::pair<JunctionKey, JunctionData>> result;
    result.reserve(batch.size());
    std::mutex result_mutex;

    // Create a function to process a single junction
    auto process_junction = [&](size_t idx) {
        auto& [key, data] = batch[idx];

        auto rseq_it = seqs.find(data.jstart_region);
        auto qseq_it = seqs.find(data.jend_region);

        if (rseq_it == seqs.end() || qseq_it == seqs.end()) {
            // Skip junctions without sequences (matches Python behavior)
            return;
        }

        const std::string& rseq = rseq_it->second;
        const std::string& qseq = qseq_it->second;

        // Each thread needs its own aligner (minimap2 structures aren't thread-safe)
        SelfAligner aligner(scoring, params);
        auto hit = aligner.align(rseq, qseq);

        if (!hit) {
            // Skip junctions without self-alignment hit (matches Python behavior)
            return;
        }

        // Apply the same filtering as get_self_aligned_introns
        int64_t intron_len = key.length();
        bool skip = false;

        if (overhang * 2 >= intron_len) {
            if (hit->r_st == intron_len) {
                skip = true;
            } else if (overhang > intron_len) {
                if (hit->r_st - hit->q_st == intron_len) {
                    skip = true;
                }
            } else if (hit->r_st >= overhang && hit->q_en <= overhang) {
                skip = true;
            }
        }

        if (skip) {
            // Skip filtered junctions (matches Python behavior)
            return;
        } else {
            JunctionData new_data = data;
            new_data.hit = hit;
            std::lock_guard<std::mutex> lock(result_mutex);
            result.emplace_back(key, std::move(new_data));
        }
    };

    if (num_threads <= 1) {
        // Sequential processing
        for (size_t i = 0; i < batch.size(); ++i) {
            process_junction(i);
        }
    } else {
        // Parallel processing using thread pool
        ThreadPool pool(static_cast<size_t>(num_threads));
        std::vector<std::future<void>> futures;
        futures.reserve(batch.size());

        for (size_t i = 0; i < batch.size(); ++i) {
            futures.push_back(pool.submit(process_junction, i));
        }

        // Wait for all tasks to complete
        for (auto& f : futures) {
            f.get();
        }
    }

    return result;
}

} // namespace eastr
