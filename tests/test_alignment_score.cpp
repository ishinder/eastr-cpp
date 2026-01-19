/**
 * Test binary to output detailed alignment information for score comparison.
 *
 * Usage: ./test_alignment_score <seq_file>
 *
 * Where seq_file contains two lines:
 *   Line 1: reference sequence (5' flanking)
 *   Line 2: query sequence (3' flanking)
 *
 * Output: Key-value pairs with alignment details
 */

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "eastr/self_aligner.hpp"
#include "eastr/types.hpp"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <seq_file>" << std::endl;
        std::cerr << "  seq_file: text file with rseq on line 1, qseq on line 2" << std::endl;
        return 1;
    }

    // Read sequences from file
    std::ifstream infile(argv[1]);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open file " << argv[1] << std::endl;
        return 1;
    }

    std::string rseq, qseq;
    std::getline(infile, rseq);
    std::getline(infile, qseq);
    infile.close();

    if (rseq.empty() || qseq.empty()) {
        std::cerr << "Error: Empty sequences" << std::endl;
        return 1;
    }

    // Use default EASTR scoring parameters
    eastr::ScoringParams scoring;
    scoring.match_score = 3;
    scoring.mismatch_penalty = 4;
    scoring.gap_open1 = 12;
    scoring.gap_extend1 = 2;
    scoring.gap_open2 = 32;
    scoring.gap_extend2 = 1;

    eastr::AlgorithmParams params;
    params.kmer_size = 3;
    params.window_size = 2;
    params.min_chain_score = 25;

    // Create aligner and align
    eastr::SelfAligner aligner(scoring, params);
    auto hit = aligner.align(rseq, qseq);

    if (!hit) {
        std::cout << "hit: none" << std::endl;
        return 0;
    }

    // Output all alignment details
    std::cout << "r_st: " << hit->r_st << std::endl;
    std::cout << "r_en: " << hit->r_en << std::endl;
    std::cout << "q_st: " << hit->q_st << std::endl;
    std::cout << "q_en: " << hit->q_en << std::endl;
    std::cout << "strand: " << hit->strand << std::endl;
    std::cout << "dp_score: " << hit->score << std::endl;
    std::cout << "mlen: " << hit->mlen << std::endl;
    std::cout << "cigar: " << hit->cigar << std::endl;
    std::cout << "is_primary: " << (hit->is_primary ? "true" : "false") << std::endl;

    // Also output the raw score calculation details
    std::cout << "scoring_match: " << scoring.match_score << std::endl;
    std::cout << "scoring_mismatch: " << scoring.mismatch_penalty << std::endl;
    std::cout << "scoring_gap_open1: " << scoring.gap_open1 << std::endl;
    std::cout << "scoring_gap_ext1: " << scoring.gap_extend1 << std::endl;

    // Calculate expected score from CIGAR by comparing sequences
    // This mimics what Python does with the cs tag
    int r_pos = hit->r_st;
    int q_pos = hit->q_st;
    int matches = 0;
    int mismatches = 0;
    int gap_penalty_python = 0;  // Python formula: gap_open + (len-1)*gap_ext
    int gap_penalty_mm2 = 0;     // mm2 formula: gap_open + len*gap_ext

    // Parse CIGAR
    std::string cigar = hit->cigar;
    size_t i = 0;
    while (i < cigar.length()) {
        size_t j = i;
        while (j < cigar.length() && std::isdigit(cigar[j])) j++;
        int len = std::stoi(cigar.substr(i, j - i));
        char op = cigar[j];

        if (op == 'M') {
            // Count matches and mismatches
            for (int k = 0; k < len; k++) {
                if (rseq[r_pos + k] == qseq[q_pos + k]) {
                    matches++;
                } else {
                    mismatches++;
                }
            }
            r_pos += len;
            q_pos += len;
        } else if (op == 'I') {
            // Insertion in query
            int py_pen = std::min(scoring.gap_open1 + (len - 1) * scoring.gap_extend1,
                                  scoring.gap_open2 + (len - 1) * scoring.gap_extend2);
            int mm2_pen = std::min(scoring.gap_open1 + len * scoring.gap_extend1,
                                   scoring.gap_open2 + len * scoring.gap_extend2);
            gap_penalty_python += py_pen;
            gap_penalty_mm2 += mm2_pen;
            q_pos += len;
        } else if (op == 'D') {
            // Deletion in query
            int py_pen = std::min(scoring.gap_open1 + (len - 1) * scoring.gap_extend1,
                                  scoring.gap_open2 + (len - 1) * scoring.gap_extend2);
            int mm2_pen = std::min(scoring.gap_open1 + len * scoring.gap_extend1,
                                   scoring.gap_open2 + len * scoring.gap_extend2);
            gap_penalty_python += py_pen;
            gap_penalty_mm2 += mm2_pen;
            r_pos += len;
        }
        i = j + 1;
    }

    int expected_score_python = matches * scoring.match_score
                               - mismatches * scoring.mismatch_penalty
                               - gap_penalty_python;
    int expected_score_mm2 = matches * scoring.match_score
                            - mismatches * scoring.mismatch_penalty
                            - gap_penalty_mm2;

    std::cout << "\n--- Score Calculation from CIGAR ---" << std::endl;
    std::cout << "matches: " << matches << std::endl;
    std::cout << "mismatches: " << mismatches << std::endl;
    std::cout << "gap_penalty_python_formula: " << gap_penalty_python << std::endl;
    std::cout << "gap_penalty_mm2_formula: " << gap_penalty_mm2 << std::endl;
    std::cout << "expected_score_python_formula: " << expected_score_python << std::endl;
    std::cout << "expected_score_mm2_formula: " << expected_score_mm2 << std::endl;
    std::cout << "actual_dp_score: " << hit->score << std::endl;
    std::cout << "diff_from_mm2_expected: " << (expected_score_mm2 - hit->score) << std::endl;

    return 0;
}
