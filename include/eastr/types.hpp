#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace eastr {

// Strand enumeration
enum class Strand : char { Plus = '+', Minus = '-', Unknown = '.' };

inline char to_char(Strand s) {
    return static_cast<char>(s);
}

inline Strand strand_from_char(char c) {
    switch (c) {
        case '+': return Strand::Plus;
        case '-': return Strand::Minus;
        default: return Strand::Unknown;
    }
}

// Scoring parameters (matching Python's minimap2/mappy parameters)
struct ScoringParams {
    int match_score = 3;        // -A
    int mismatch_penalty = 4;   // -B
    int gap_open1 = 12;         // -O[0]
    int gap_extend1 = 2;        // -E[0]
    int gap_open2 = 32;         // -O[1]
    int gap_extend2 = 1;        // -E[1]
    int ambiguous_score = 1;    // --scoreN

    // Convert to array format for compatibility
    std::array<int, 7> to_array() const {
        return {match_score, mismatch_penalty, gap_open1, gap_extend1,
                gap_open2, gap_extend2, ambiguous_score};
    }
};

// Algorithm parameters
struct AlgorithmParams {
    int overhang = 50;                    // -o: length of flanking sequence
    int anchor = 7;                       // -a: minimum anchor length
    int min_duplicate_exon_length = 27;   // --min_duplicate_exon_length
    int min_junc_score = 1;               // --min_junc_score
    int bt2_k = 10;                       // --bt2_k: min distinct alignments
    int kmer_size = 3;                    // -k: k-mer length for alignment
    int window_size = 2;                  // -w: minimizer window size
    int min_chain_score = 25;             // -m: minimum chain score
};

// Configuration for the entire pipeline
struct Config {
    // Input paths (mutually exclusive)
    std::string gtf_path;
    std::string bed_path;      // Single file or list file
    std::string bam_path;      // Single file or list file

    // Required paths
    std::string reference_fasta;
    std::string bowtie2_index;

    // Parameters
    ScoringParams scoring;
    AlgorithmParams params;

    // Output options
    std::string out_original_junctions;
    std::string out_removed_junctions = "stdout";
    std::string out_kept_junctions;  // Non-spurious junctions
    std::string out_filtered_bam;
    std::string filtered_bam_suffix = "_EASTR_filtered";

    // Processing options
    int num_threads = 1;
    int batch_size = 20000;  // Number of junctions to batch for streaming
    int prefetch_buffer_size = 10000;  // Prefetch buffer size for BAM reading
    bool verbose = false;
    bool write_removed_alignments = false;

    // Trusted junctions (optional)
    std::string trusted_bed;
};

} // namespace eastr
