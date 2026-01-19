#pragma once

#include "eastr/junction.hpp"
#include "eastr/multi_bam_reader.hpp"

#include <htslib/sam.h>

#include <map>
#include <string>
#include <vector>
#include <utility>

namespace eastr {

/**
 * Accumulates junction scores from streaming alignments.
 *
 * Tracks "active" junctions (those at recent positions) and detects
 * when a junction is "complete" (all files have advanced past its position).
 *
 * Memory usage: O(W) where W = active window size (junctions at recent positions)
 */
class JunctionAccumulator {
public:
    /**
     * Create accumulator for the given number of BAM files.
     * @param num_files Number of BAM files being processed
     * @param sample_names Names for each sample (for JunctionData)
     */
    JunctionAccumulator(size_t num_files, const std::vector<std::string>& sample_names);

    /**
     * Process an alignment record, extracting and accumulating any junctions.
     * @param rec Alignment record from MultiBamReader
     * @param header BAM header (for chromosome name lookup)
     */
    void process_alignment(const MultiBamReader::AlignmentRecord& rec, sam_hdr_t* header);

    /**
     * Get junctions that are now complete (all files have advanced past them).
     * @return Vector of completed junctions (moved out, no longer tracked)
     */
    std::vector<std::pair<JunctionKey, JunctionData>> get_completed_junctions();

    /**
     * Flush all remaining junctions for the current chromosome.
     * Call this when chromosome changes or at end of processing.
     * @return All remaining active junctions
     */
    std::vector<std::pair<JunctionKey, JunctionData>> flush();

    /**
     * Notify that a file has reached a new chromosome.
     * @param file_idx Which file changed chromosomes
     * @param new_tid New chromosome ID (-1 if file exhausted)
     */
    void notify_chromosome_change(int file_idx, int new_tid);

    /**
     * Get count of currently active junctions.
     */
    size_t active_junction_count() const { return active_junctions_.size(); }

    /**
     * Get total junctions seen so far.
     */
    size_t total_junctions_seen() const { return total_junctions_; }

private:
    /**
     * Extract junctions from a spliced alignment's CIGAR string.
     * @param rec Alignment record
     * @param header BAM header
     * @return Vector of (junction_key, score) pairs
     */
    std::vector<std::pair<JunctionKey, int>> extract_junctions_from_alignment(
        const MultiBamReader::AlignmentRecord& rec,
        sam_hdr_t* header);

    /**
     * Check if a position is complete (all files have advanced past it).
     */
    bool is_position_complete(const std::string& chrom, int64_t end_pos) const;

    /**
     * Update minimum position for a file.
     */
    void update_file_position(int file_idx, const std::string& chrom, int64_t pos);

    // Per-file tracking: current chromosome and minimum position
    struct FileState {
        std::string current_chrom;
        int64_t min_pos = -1;
        bool exhausted = false;
    };
    std::vector<FileState> file_states_;

    // Sample names for creating JunctionData
    std::vector<std::string> sample_names_;

    // Active junctions (keyed by JunctionKey for deduplication)
    // Using std::map for ordered iteration (can emit in coordinate order)
    std::map<JunctionKey, JunctionData> active_junctions_;

    // Statistics
    size_t total_junctions_ = 0;

    // Current chromosome being processed (for completion detection)
    std::string current_chrom_;
};

} // namespace eastr
