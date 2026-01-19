#pragma once

#include "junction.hpp"
#include "types.hpp"

#include <memory>
#include <optional>
#include <string>

namespace eastr {

class SelfAligner {
public:
    explicit SelfAligner(const ScoringParams& scoring, const AlgorithmParams& params);
    ~SelfAligner();

    // Non-copyable, movable
    SelfAligner(const SelfAligner&) = delete;
    SelfAligner& operator=(const SelfAligner&) = delete;
    SelfAligner(SelfAligner&&) noexcept;
    SelfAligner& operator=(SelfAligner&&) noexcept;

    // Align query sequence to reference sequence
    // Returns nullopt if no significant alignment found
    std::optional<AlignmentHit> align(const std::string& ref_seq,
                                       const std::string& query_seq);

    // Calculate alignment score from hit
    int calc_alignment_score(const AlignmentHit& hit) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl_;
};

// Process all junctions to find self-alignments
// Returns junctions that have significant self-alignments
// Modifies input junctions to add hit information
JunctionMap get_self_aligned_introns(
    const JunctionMap& introns,
    const SequenceCache& seqs,
    int overhang,
    const ScoringParams& scoring,
    const AlgorithmParams& params);

// Process batch of junctions in parallel to find self-alignments
// Modifies junctions in-place to add hit information
// Returns vector with all junctions (hit information added where applicable)
std::vector<std::pair<JunctionKey, JunctionData>> get_self_aligned_introns_batch(
    std::vector<std::pair<JunctionKey, JunctionData>>& batch,
    const SequenceCache& seqs,
    int overhang,
    const ScoringParams& scoring,
    const AlgorithmParams& params,
    int num_threads = 1);

} // namespace eastr
