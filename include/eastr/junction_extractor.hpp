#pragma once

#include "junction.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace eastr {

class JunctionExtractor {
public:
    explicit JunctionExtractor(const std::string& bam_path);
    ~JunctionExtractor();

    // Non-copyable, movable
    JunctionExtractor(const JunctionExtractor&) = delete;
    JunctionExtractor& operator=(const JunctionExtractor&) = delete;
    JunctionExtractor(JunctionExtractor&&) noexcept;
    JunctionExtractor& operator=(JunctionExtractor&&) noexcept;

    // Extract all junctions from the BAM file
    // Returns map of junction key -> (sample_name, score)
    std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash>
    extract();

    // Write junctions to BED file
    void write_bed(const std::string& output_path);

    // Get the sample name (derived from BAM filename)
    const std::string& sample_name() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl_;
};

// Extract junctions from multiple BAM files in parallel
JunctionMap extract_junctions_multi_bam(
    const std::vector<std::string>& bam_paths,
    const std::vector<std::string>& output_paths,
    int num_threads);

} // namespace eastr
