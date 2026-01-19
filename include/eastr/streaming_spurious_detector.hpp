#pragma once

#include "eastr/junction.hpp"
#include "eastr/fasta_index.hpp"
#include "eastr/types.hpp"

#include <memory>
#include <string>
#include <vector>
#include <utility>

namespace eastr {

/**
 * Streaming spurious junction detector.
 *
 * Processes junctions in batches as they complete from the accumulator,
 * rather than loading all junctions into memory at once.
 *
 * Memory usage: O(batch_size) for junctions being tested
 */
class StreamingSpuriousDetector {
public:
    /**
     * Create a streaming spurious detector.
     * @param fasta Reference FASTA index for sequence retrieval
     * @param bt2_index Path to bowtie2 index
     * @param scoring Scoring parameters for alignment
     * @param params Algorithm parameters
     * @param batch_size Number of junctions to batch before testing
     * @param num_threads Number of threads for parallel processing
     */
    StreamingSpuriousDetector(
        FastaIndex& fasta,
        const std::string& bt2_index,
        const ScoringParams& scoring,
        const AlgorithmParams& params,
        int batch_size = 20000,
        int num_threads = 1);

    ~StreamingSpuriousDetector();

    // Non-copyable
    StreamingSpuriousDetector(const StreamingSpuriousDetector&) = delete;
    StreamingSpuriousDetector& operator=(const StreamingSpuriousDetector&) = delete;

    /**
     * Add a completed junction to the testing queue.
     * @param key Junction key
     * @param data Junction data with accumulated scores
     */
    void add_junction(JunctionKey key, JunctionData data);

    /**
     * Check if the batch is full and ready for processing.
     */
    bool batch_full() const { return batch_.size() >= static_cast<size_t>(batch_size_); }

    /**
     * Get current batch size.
     */
    size_t current_batch_size() const { return batch_.size(); }

    /**
     * Process the current batch and return spurious junctions.
     * Clears the batch after processing.
     * @param is_bam_input True if processing BAM input (affects score thresholds)
     * @return Vector of spurious junctions from this batch
     */
    std::vector<std::pair<JunctionKey, JunctionData>> process_batch(bool is_bam_input);

    /**
     * Flush any remaining junctions in the batch.
     * @param is_bam_input True if processing BAM input
     * @return Vector of spurious junctions from final batch
     */
    std::vector<std::pair<JunctionKey, JunctionData>> flush(bool is_bam_input);

    /**
     * Get statistics.
     */
    size_t total_junctions_processed() const { return total_processed_; }
    size_t total_spurious_found() const { return total_spurious_; }

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl_;

    // Current batch of junctions awaiting processing
    std::vector<std::pair<JunctionKey, JunctionData>> batch_;

    int batch_size_;
    int num_threads_;

    // Statistics
    size_t total_processed_ = 0;
    size_t total_spurious_ = 0;
};

} // namespace eastr
