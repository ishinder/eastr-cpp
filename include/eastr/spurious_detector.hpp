#pragma once

#include "junction.hpp"
#include "types.hpp"
#include "fasta_index.hpp"
#include "self_aligner.hpp"
#include "bowtie2_runner.hpp"

#include <memory>
#include <optional>
#include <string>
#include <unordered_set>

namespace eastr {

class SpuriousDetector {
public:
    SpuriousDetector(const ScoringParams& scoring,
                     const AlgorithmParams& params,
                     const std::string& bt2_index,
                     const std::string& ref_fasta,
                     int num_threads);

    ~SpuriousDetector();

    // Main entry point: identify spurious junctions
    JunctionMap detect_spurious(
        const JunctionMap& junctions,
        bool is_bam_input,
        const std::optional<std::string>& trusted_bed_path = std::nullopt);

private:
    ScoringParams scoring_;
    AlgorithmParams params_;
    std::string bt2_index_;
    std::unique_ptr<FastaIndex> fasta_;
    int num_threads_;

    // Two-anchor alignment check
    static bool is_two_anchor_alignment(
        const AlignmentHit& hit,
        int overhang,
        int anchor);

    // Check if junction is spurious based on alignment patterns
    bool is_spurious_alignment(
        const JunctionKey& key,
        const JunctionData& data,
        const SequenceCache& seqs) const;

    // Linear distance between two equal-length strings
    static int linear_distance(const std::string& s1, const std::string& s2);

    // Load trusted junctions from BED file
    std::unordered_set<JunctionKey, JunctionKeyHash>
    load_trusted_junctions(const std::string& bed_path);
};

} // namespace eastr
