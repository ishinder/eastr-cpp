#pragma once

#include "eastr/junction.hpp"

#include <string>
#include <vector>

namespace eastr {

/**
 * Extract junctions from multiple BAM files in parallel.
 *
 * Supports two levels of parallelism:
 * 1. Inter-file: Multiple BAM files processed concurrently
 * 2. Intra-file: Each BAM file processed by genomic region (tiled)
 *
 * This avoids the bottleneck of sequential streaming through merged alignment
 * records and enables full CPU utilization even with few BAM files.
 */
class ParallelBamExtractor {
public:
    struct Options {
        int num_threads = 1;
        bool verbose = false;
        int64_t tile_size = 10000000;  // 10MB tiles for intra-file parallelism
        int htslib_threads = 2;         // Threads for BGZF decompression per file
    };

    /**
     * Extract junctions from multiple BAM files in parallel.
     *
     * Uses both inter-file and intra-file parallelism:
     * - If num_threads >= num_files: one thread per file + htslib decompression
     * - If num_threads > num_files: uses tiled processing within files
     *
     * @param bam_paths Paths to BAM files
     * @param num_threads Total number of threads to use
     * @param verbose Print progress information
     * @return Merged junction map with scores from all samples
     */
    static JunctionMap extract_junctions_parallel(
        const std::vector<std::string>& bam_paths,
        int num_threads,
        bool verbose = false);

    /**
     * Extract junctions from a single BAM file (sequential).
     *
     * @param bam_path Path to BAM file
     * @param sample_name Name to use for this sample in junction data
     * @param htslib_threads Threads for BGZF decompression (default: 2)
     * @return Junction map for this file
     */
    static JunctionMap extract_junctions_from_bam(
        const std::string& bam_path,
        const std::string& sample_name,
        int htslib_threads = 2);

    /**
     * Extract junctions from a single BAM file using tiled parallelism.
     *
     * Divides the BAM file into genomic regions and processes them in parallel.
     * Requires BAM index (.bai) file.
     *
     * @param bam_path Path to BAM file
     * @param sample_name Name to use for this sample in junction data
     * @param num_threads Number of threads for parallel region processing
     * @param tile_size Size of each genomic tile in base pairs
     * @param verbose Print progress information
     * @return Junction map for this file
     */
    static JunctionMap extract_junctions_from_bam_tiled(
        const std::string& bam_path,
        const std::string& sample_name,
        int num_threads,
        int64_t tile_size = 10000000,
        bool verbose = false);

    /**
     * Extract junctions from a specific genomic region.
     *
     * @param bam_path Path to BAM file
     * @param sample_name Name to use for this sample in junction data
     * @param region Region string (e.g., "chr1:1000000-2000000")
     * @return Junction map for this region
     */
    static JunctionMap extract_junctions_from_region(
        const std::string& bam_path,
        const std::string& sample_name,
        const std::string& region);

    /**
     * Merge multiple junction maps into one.
     * Combines scores and sample information.
     *
     * @param maps Vector of junction maps to merge
     * @return Merged junction map
     */
    static JunctionMap merge_junction_maps(std::vector<JunctionMap>& maps);

private:
    /**
     * Check if BAM index exists for given BAM file.
     */
    static bool has_bam_index(const std::string& bam_path);
};

} // namespace eastr
