#pragma once

#include "types.hpp"

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace eastr {

// Junction key: (chrom, start, end, strand)
// Using 0-based, half-open coordinates internally (BED-style)
struct JunctionKey {
    std::string chrom;
    int64_t start;      // 0-based, inclusive
    int64_t end;        // 0-based, exclusive
    Strand strand;

    bool operator==(const JunctionKey& other) const {
        return chrom == other.chrom && start == other.start &&
               end == other.end && strand == other.strand;
    }

    bool operator<(const JunctionKey& other) const {
        if (chrom != other.chrom) return chrom < other.chrom;
        if (start != other.start) return start < other.start;
        if (end != other.end) return end < other.end;
        return strand < other.strand;
    }

    // Intron length
    int64_t length() const { return end - start; }

    // Convert to string for display/debugging
    std::string to_string() const;
};

// Custom hash for JunctionKey
struct JunctionKeyHash {
    size_t operator()(const JunctionKey& key) const {
        size_t h1 = std::hash<std::string>{}(key.chrom);
        size_t h2 = std::hash<int64_t>{}(key.start);
        size_t h3 = std::hash<int64_t>{}(key.end);
        size_t h4 = std::hash<char>{}(static_cast<char>(key.strand));
        // Combine hashes using boost-style hash_combine
        size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

// Sample information for a junction
struct SampleInfo {
    std::string sample_name;
    int score;  // Number of supporting reads

    bool operator==(const SampleInfo& other) const {
        return sample_name == other.sample_name && score == other.score;
    }
};

// Alignment hit from self-alignment (ksw2)
struct AlignmentHit {
    int r_st = 0;       // Reference (5' flanking) start
    int r_en = 0;       // Reference (5' flanking) end
    int q_st = 0;       // Query (3' flanking) start
    int q_en = 0;       // Query (3' flanking) end
    int strand = 1;     // 1 for forward, -1 for reverse
    int score = 0;      // Alignment score
    std::string cigar;  // CIGAR string
    bool is_primary = true;

    // Match length (number of matching positions)
    int mlen = 0;
};

// Complete junction data
struct JunctionData {
    int total_score = 0;                    // Sum of scores across samples
    std::vector<SampleInfo> samples;        // Per-sample information

    // Flanking region identifiers (format: "chr:start-end")
    std::string jstart_region;              // 5' flank region
    std::string jend_region;                // 3' flank region

    // Self-alignment results (populated during analysis)
    std::optional<AlignmentHit> hit;

    // Bowtie2 alignment counts
    int seq1_count = 0;     // Alignments for ref overlap region
    int seq2_count = 0;     // Alignments for query overlap region
    int seqh_count = 0;     // Alignments for hybrid sequence

    // GTF-specific fields
    std::string gene_id;
    std::vector<std::string> transcript_ids;

    // Add a sample's contribution (consolidates by sample name)
    void add_sample(const std::string& name, int score) {
        // Check if sample already exists and consolidate
        for (auto& sample : samples) {
            if (sample.sample_name == name) {
                sample.score += score;
                total_score += score;
                return;
            }
        }
        // New sample
        samples.push_back({name, score});
        total_score += score;
    }

    // Merge another junction's data into this one
    void merge(const JunctionData& other) {
        for (const auto& sample : other.samples) {
            samples.push_back(sample);
        }
        total_score += other.total_score;
    }
};

// Main junction container type
using JunctionMap = std::unordered_map<JunctionKey, JunctionData, JunctionKeyHash>;

// Flanking sequence cache: region string -> sequence
using SequenceCache = std::unordered_map<std::string, std::string>;

// Helper function to create region string
inline std::string make_region_string(const std::string& chrom, int64_t start, int64_t end) {
    return chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
}

// Parse region string back to components
bool parse_region_string(const std::string& region, std::string& chrom,
                         int64_t& start, int64_t& end);

} // namespace eastr
