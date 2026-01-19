#pragma once

#include <htslib/sam.h>

#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace eastr {

/**
 * Coordinate-sorted multi-BAM reader.
 *
 * Iterates multiple BAM files simultaneously in coordinate order,
 * like tiebrush/stringtie's approach. This enables streaming processing
 * without loading all data into memory.
 *
 * Supports per-file prefetching for improved I/O parallelism.
 *
 * Memory usage: O(N * buffer_size) where N = number of BAM files
 */
class MultiBamReader {
public:
    /**
     * Alignment record with source file information.
     */
    struct AlignmentRecord {
        bam1_t* aln;            // The alignment (owned by MultiBamReader)
        int file_index;         // Which file this came from
        std::string sample_name; // Sample name (derived from filename)
        int tid;                // Chromosome/contig ID
        int64_t pos;            // 0-based position
    };

    /**
     * Create a multi-BAM reader.
     * @param bam_paths Paths to BAM files (must be coordinate-sorted)
     * @throws std::runtime_error if any file is not coordinate-sorted
     */
    explicit MultiBamReader(const std::vector<std::string>& bam_paths);

    ~MultiBamReader();

    // Non-copyable
    MultiBamReader(const MultiBamReader&) = delete;
    MultiBamReader& operator=(const MultiBamReader&) = delete;

    /**
     * Get the next alignment in coordinate order across all files.
     * @return Pointer to alignment record, or nullptr if all files exhausted
     */
    AlignmentRecord* next();

    /**
     * Get the shared BAM header.
     * All BAM files should have compatible headers.
     */
    sam_hdr_t* header() const { return header_; }

    /**
     * Get number of input files.
     */
    size_t num_files() const { return files_.size(); }

    /**
     * Enable per-file prefetching with background threads.
     * @param buffer_size Number of alignments to buffer per file (default: 10000)
     */
    void enable_prefetching(size_t buffer_size = 10000);

    /**
     * Disable prefetching and stop all background threads.
     */
    void disable_prefetching();

    /**
     * Check if prefetching is enabled.
     */
    bool is_prefetching_enabled() const { return prefetching_enabled_; }

    /**
     * Get current total buffer occupancy (for diagnostics).
     */
    size_t buffer_occupancy() const;

private:
    /**
     * Validate that a BAM file is coordinate-sorted.
     */
    void validate_sorted(samFile* fp, sam_hdr_t* hdr, const std::string& path);

    /**
     * Read the next record from a specific file (synchronous).
     */
    bool read_next(int file_idx);

    /**
     * Synchronous next() implementation (original logic).
     */
    AlignmentRecord* next_sync();

    /**
     * Prefetched next() using per-file buffers.
     */
    AlignmentRecord* next_prefetch();

    /**
     * Entry in the priority queue for coordinate-sorted merging.
     */
    struct QueueEntry {
        int file_idx;
        int tid;
        int64_t pos;

        // Min-heap: smallest (tid, pos) first
        bool operator>(const QueueEntry& other) const {
            if (tid != other.tid) return tid > other.tid;
            if (pos != other.pos) return pos > other.pos;
            return file_idx > other.file_idx;
        }
    };

    std::vector<samFile*> files_;
    std::vector<bam1_t*> records_;       // One buffered record per file (sync mode)
    std::vector<std::string> sample_names_;
    std::vector<bool> exhausted_;

    sam_hdr_t* header_ = nullptr;        // Shared header

    std::priority_queue<QueueEntry, std::vector<QueueEntry>,
                        std::greater<QueueEntry>> queue_;

    AlignmentRecord current_record_;     // Returned by next()
    bam1_t* sync_current_aln_ = nullptr;      // Owned copy for sync mode
    bam1_t* prefetch_current_aln_ = nullptr;  // Owned copy for prefetch mode

    // ========== Per-file prefetch state ==========

    bool prefetching_enabled_ = false;
    size_t per_file_buffer_size_ = 1000;  // Buffer size per file

    // Buffered record for prefetching (includes deep copy of alignment data)
    struct BufferedRecord {
        bam1_t* aln;           // Deep copy of alignment
        int file_index;
        std::string sample_name;
        int tid;
        int64_t pos;
        bool valid;

        BufferedRecord() : aln(nullptr), file_index(-1), tid(-1), pos(-1), valid(false) {}
    };

    // Per-file buffer and thread state
    struct PerFileState {
        std::vector<BufferedRecord> buffer;
        std::atomic<size_t> read_pos{0};
        std::atomic<size_t> write_pos{0};
        std::mutex mutex;
        std::condition_variable not_empty;
        std::condition_variable not_full;
        std::thread thread;
        std::atomic<bool> eof{false};
        std::atomic<bool> stop{false};
        size_t buffer_size = 0;

        PerFileState() = default;
        PerFileState(const PerFileState&) = delete;
        PerFileState& operator=(const PerFileState&) = delete;
    };

    std::vector<std::unique_ptr<PerFileState>> per_file_states_;

    // Current record from each file's buffer (for priority queue)
    std::vector<BufferedRecord> current_from_file_;

    /**
     * Per-file prefetch worker thread function.
     */
    void per_file_prefetch_worker(int file_idx);

    /**
     * Pop next record from a file's buffer.
     * Returns false if buffer is empty and EOF reached.
     */
    bool pop_from_file_buffer(int file_idx, BufferedRecord& out);

    /**
     * Peek at the next record from a file's buffer (don't consume).
     * Returns false if buffer is empty and EOF reached.
     */
    bool peek_file_buffer(int file_idx, int& tid, int64_t& pos);

    /**
     * Refill the priority queue entry for a specific file.
     */
    bool refill_queue_entry(int file_idx);
};

} // namespace eastr
