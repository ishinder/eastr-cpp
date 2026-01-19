#include "eastr/fasta_index.hpp"

#include <htslib/faidx.h>

#include <algorithm>
#include <cstring>
#include <future>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_set>
#include <vector>

namespace eastr {

struct FastaIndex::Impl {
    faidx_t* fai = nullptr;
    std::string path;

    ~Impl() {
        if (fai) {
            fai_destroy(fai);
        }
    }
};

FastaIndex::FastaIndex(const std::string& fasta_path)
    : pImpl_(std::make_unique<Impl>()) {
    pImpl_->path = fasta_path;

    // Ensure index exists
    ensure_indexed(fasta_path);

    // Load the index
    pImpl_->fai = fai_load(fasta_path.c_str());
    if (!pImpl_->fai) {
        throw std::runtime_error("Failed to load FASTA index for: " + fasta_path);
    }
}

FastaIndex::~FastaIndex() = default;

FastaIndex::FastaIndex(FastaIndex&&) noexcept = default;
FastaIndex& FastaIndex::operator=(FastaIndex&&) noexcept = default;

std::string FastaIndex::fetch(const std::string& chrom, int64_t start, int64_t end) {
    if (!pImpl_->fai) {
        throw std::runtime_error("FASTA index not loaded");
    }

    // faidx uses 0-based coordinates internally
    // Input is 1-based inclusive (samtools convention), convert to 0-based
    int64_t start_0based = start - 1;
    int64_t end_0based = end - 1;

    hts_pos_t len;
    char* seq = faidx_fetch_seq64(pImpl_->fai, chrom.c_str(), start_0based, end_0based, &len);

    if (!seq) {
        throw std::runtime_error("Failed to fetch sequence for region: " +
                                 chrom + ":" + std::to_string(start) + "-" + std::to_string(end));
    }

    std::string result(seq, len);
    free(seq);

    // Convert to uppercase
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);

    return result;
}

std::string FastaIndex::fetch_region(const std::string& region) {
    std::string chrom;
    int64_t start, end;

    if (!parse_region_string(region, chrom, start, end)) {
        throw std::runtime_error("Invalid region string: " + region);
    }

    return fetch(chrom, start, end);
}

int64_t FastaIndex::get_length(const std::string& chrom) const {
    if (!pImpl_->fai) {
        throw std::runtime_error("FASTA index not loaded");
    }

    int len = faidx_seq_len(pImpl_->fai, chrom.c_str());
    if (len < 0) {
        throw std::runtime_error("Chromosome not found in FASTA: " + chrom);
    }

    return static_cast<int64_t>(len);
}

std::unordered_map<std::string, int64_t> FastaIndex::get_chrom_sizes() const {
    if (!pImpl_->fai) {
        throw std::runtime_error("FASTA index not loaded");
    }

    std::unordered_map<std::string, int64_t> sizes;
    int nseq = faidx_nseq(pImpl_->fai);

    for (int i = 0; i < nseq; ++i) {
        const char* name = faidx_iseq(pImpl_->fai, i);
        int len = faidx_seq_len(pImpl_->fai, name);
        sizes[name] = static_cast<int64_t>(len);
    }

    return sizes;
}

void FastaIndex::ensure_indexed(const std::string& fasta_path) {
    std::string index_path = fasta_path + ".fai";

    // Check if index exists
    FILE* f = fopen(index_path.c_str(), "r");
    if (f) {
        fclose(f);
        return;  // Index exists
    }

    // Create index
    if (fai_build(fasta_path.c_str()) != 0) {
        throw std::runtime_error("Failed to create FASTA index for: " + fasta_path);
    }
}

const std::string& FastaIndex::get_path() const {
    return pImpl_->path;
}

// Region to fetch with metadata for sorting
struct RegionToFetch {
    std::string chrom;
    int64_t start;
    int64_t end;
    std::string region_string;

    // Sort by chromosome name, then by start position
    bool operator<(const RegionToFetch& other) const {
        if (chrom != other.chrom) return chrom < other.chrom;
        return start < other.start;
    }
};

SequenceCache get_flanking_subsequences(
    JunctionMap& junctions,
    const std::unordered_map<std::string, int64_t>& chrom_sizes,
    int overhang,
    FastaIndex& fasta) {

    SequenceCache seqs;
    std::vector<JunctionKey> to_remove;
    std::vector<RegionToFetch> regions_to_fetch;

    // First pass: collect all regions to fetch and update junction data
    for (auto& [key, data] : junctions) {
        auto chrom_it = chrom_sizes.find(key.chrom);
        if (chrom_it == chrom_sizes.end()) {
            // Chromosome not in FASTA
            to_remove.push_back(key);
            continue;
        }

        int64_t max_length = chrom_it->second;

        // Calculate flanking regions (1-based coordinates for samtools)
        int64_t rstart = std::max(key.start - overhang + 1, int64_t(1));
        int64_t rend = std::min(key.start + overhang, max_length);
        int64_t qstart = std::max(key.end - overhang + 1, int64_t(1));
        int64_t qend = std::min(key.end + overhang, max_length);

        // Check bounds
        if (rstart > max_length || qstart > max_length) {
            to_remove.push_back(key);
            continue;
        }

        // Create region strings
        std::string r1 = make_region_string(key.chrom, rstart, rend);
        std::string r2 = make_region_string(key.chrom, qstart, qend);

        data.jstart_region = r1;
        data.jend_region = r2;

        // Collect regions to fetch (will deduplicate via the cache)
        regions_to_fetch.push_back({key.chrom, rstart, rend, r1});
        regions_to_fetch.push_back({key.chrom, qstart, qend, r2});
    }

    // Remove invalid junctions
    for (const auto& key : to_remove) {
        junctions.erase(key);
    }

    // Sort regions by chromosome and position for better cache locality
    std::sort(regions_to_fetch.begin(), regions_to_fetch.end());

    // Fetch sequences in sorted order (deduplicated via cache)
    for (const auto& region : regions_to_fetch) {
        if (seqs.find(region.region_string) == seqs.end()) {
            seqs[region.region_string] = fasta.fetch(region.chrom, region.start, region.end);
        }
    }

    return seqs;
}

SequenceCache get_flanking_subsequences_parallel(
    JunctionMap& junctions,
    const std::unordered_map<std::string, int64_t>& chrom_sizes,
    int overhang,
    const std::string& fasta_path,
    int num_threads) {

    // If single-threaded, use the sequential version
    if (num_threads <= 1) {
        FastaIndex fasta(fasta_path);
        return get_flanking_subsequences(junctions, chrom_sizes, overhang, fasta);
    }

    std::vector<JunctionKey> to_remove;
    std::vector<RegionToFetch> regions_to_fetch;

    // First pass: collect all regions to fetch and update junction data
    for (auto& [key, data] : junctions) {
        auto chrom_it = chrom_sizes.find(key.chrom);
        if (chrom_it == chrom_sizes.end()) {
            to_remove.push_back(key);
            continue;
        }

        int64_t max_length = chrom_it->second;

        int64_t rstart = std::max(key.start - overhang + 1, int64_t(1));
        int64_t rend = std::min(key.start + overhang, max_length);
        int64_t qstart = std::max(key.end - overhang + 1, int64_t(1));
        int64_t qend = std::min(key.end + overhang, max_length);

        if (rstart > max_length || qstart > max_length) {
            to_remove.push_back(key);
            continue;
        }

        std::string r1 = make_region_string(key.chrom, rstart, rend);
        std::string r2 = make_region_string(key.chrom, qstart, qend);

        data.jstart_region = r1;
        data.jend_region = r2;

        regions_to_fetch.push_back({key.chrom, rstart, rend, r1});
        regions_to_fetch.push_back({key.chrom, qstart, qend, r2});
    }

    // Remove invalid junctions
    for (const auto& key : to_remove) {
        junctions.erase(key);
    }

    // Sort regions by chromosome and position for better cache locality
    std::sort(regions_to_fetch.begin(), regions_to_fetch.end());

    // Deduplicate regions (keep unique region strings)
    std::vector<RegionToFetch> unique_regions;
    {
        std::unordered_set<std::string> seen;
        for (const auto& region : regions_to_fetch) {
            if (seen.insert(region.region_string).second) {
                unique_regions.push_back(region);
            }
        }
    }

    // Partition regions across threads
    size_t total_regions = unique_regions.size();
    size_t regions_per_thread = (total_regions + num_threads - 1) / num_threads;

    // Launch parallel fetches
    std::vector<std::future<SequenceCache>> futures;
    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * regions_per_thread;
        size_t end = std::min(start + regions_per_thread, total_regions);

        if (start >= end) break;

        futures.push_back(std::async(std::launch::async, [&, start, end]() {
            SequenceCache local_cache;
            FastaIndex local_fasta(fasta_path);

            for (size_t i = start; i < end; ++i) {
                const auto& region = unique_regions[i];
                local_cache[region.region_string] = local_fasta.fetch(
                    region.chrom, region.start, region.end);
            }

            return local_cache;
        }));
    }

    // Merge results
    SequenceCache seqs;
    for (auto& future : futures) {
        auto local_cache = future.get();
        for (auto& [key, value] : local_cache) {
            seqs[key] = std::move(value);
        }
    }

    return seqs;
}

} // namespace eastr
