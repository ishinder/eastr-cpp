#pragma once

#include "junction.hpp"

#include <memory>
#include <string>
#include <unordered_map>

namespace eastr {

class FastaIndex {
public:
    explicit FastaIndex(const std::string& fasta_path);
    ~FastaIndex();

    // Non-copyable, movable
    FastaIndex(const FastaIndex&) = delete;
    FastaIndex& operator=(const FastaIndex&) = delete;
    FastaIndex(FastaIndex&&) noexcept;
    FastaIndex& operator=(FastaIndex&&) noexcept;

    // Get sequence for region (1-based, inclusive coordinates for samtools compatibility)
    std::string fetch(const std::string& chrom, int64_t start, int64_t end);

    // Get sequence by region string "chr:start-end"
    std::string fetch_region(const std::string& region);

    // Get chromosome length
    int64_t get_length(const std::string& chrom) const;

    // Get all chromosome names and lengths
    std::unordered_map<std::string, int64_t> get_chrom_sizes() const;

    // Ensure index exists (create if needed)
    static void ensure_indexed(const std::string& fasta_path);

    // Get path to FASTA file (for creating per-thread instances)
    const std::string& get_path() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl_;
};

// Get flanking sequences for all junctions (batched for efficiency)
// Updates junction data with jstart_region and jend_region
// Returns cache of region -> sequence
SequenceCache get_flanking_subsequences(
    JunctionMap& junctions,
    const std::unordered_map<std::string, int64_t>& chrom_sizes,
    int overhang,
    FastaIndex& fasta);

// Parallel version of get_flanking_subsequences using multiple threads
// Each thread creates its own FastaIndex instance for thread safety
SequenceCache get_flanking_subsequences_parallel(
    JunctionMap& junctions,
    const std::unordered_map<std::string, int64_t>& chrom_sizes,
    int overhang,
    const std::string& fasta_path,
    int num_threads);

} // namespace eastr
