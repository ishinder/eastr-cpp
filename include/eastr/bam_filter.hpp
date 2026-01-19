#pragma once

#include "junction.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace eastr {

class BamFilter {
public:
    struct Options {
        bool remove_mate = true;
        bool verbose = false;
        bool write_removed = false;
        std::string removed_output_path;
        int io_threads = 1;  // Threads per BAM file for htslib I/O (BGZF compression)
    };

    explicit BamFilter(const Options& options);
    ~BamFilter();

    // Filter BAM file, removing alignments at spurious junctions
    struct FilterStats {
        int64_t total_alignments = 0;
        int64_t total_spliced = 0;
        int64_t removed_alignments = 0;
    };

    FilterStats filter(
        const std::string& input_bam,
        const std::string& spurious_bed,
        const std::string& output_bam);

private:
    Options options_;

    struct Impl;
    std::unique_ptr<Impl> pImpl_;

    // Load spurious junctions from BED
    std::unordered_set<JunctionKey, JunctionKeyHash>
    load_spurious_junctions(const std::string& bed_path);
};

// Filter multiple BAM files in parallel using a single spurious junctions BED
void filter_multi_bam(
    const std::vector<std::string>& input_bams,
    const std::string& spurious_bed,
    const std::vector<std::string>& output_bams,
    int num_threads,
    const BamFilter::Options& options);

} // namespace eastr
