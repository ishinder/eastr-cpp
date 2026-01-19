#include "eastr/parallel_bam_extractor.hpp"
#include "eastr/thread_pool.hpp"

#include <htslib/sam.h>
#include <htslib/hts.h>

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <future>
#include <iostream>
#include <mutex>
#include <sstream>

namespace fs = std::filesystem;

namespace eastr {

// Helper to extract junctions from an alignment record
static void extract_junctions_from_alignment(
    bam1_t* aln,
    sam_hdr_t* header,
    const std::string& sample_name,
    JunctionMap& junctions) {

    // Skip unmapped reads
    if (aln->core.flag & BAM_FUNMAP) {
        return;
    }

    // Get chromosome name
    const char* chrom = sam_hdr_tid2name(header, aln->core.tid);
    if (!chrom) return;

    // Check for splice junctions (N operations in CIGAR)
    uint32_t* cigar = bam_get_cigar(aln);
    bool has_splice = false;
    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
            has_splice = true;
            break;
        }
    }

    if (!has_splice) {
        return;
    }

    // Get strand from XS tag (set by aligners like HISAT2)
    Strand strand = Strand::Unknown;
    uint8_t* xs_tag = bam_aux_get(aln, "XS");
    if (xs_tag) {
        char xs_val = bam_aux2A(xs_tag);
        if (xs_val == '+') strand = Strand::Plus;
        else if (xs_val == '-') strand = Strand::Minus;
    }

    // Parse CIGAR to extract junctions
    std::string chrom_str(chrom);
    int64_t ref_pos = aln->core.pos;  // 0-based

    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CREF_SKIP) {
            // This is a splice junction (intron)
            JunctionKey key;
            key.chrom = chrom_str;
            key.start = ref_pos;
            key.end = ref_pos + len;
            key.strand = strand;

            // Add to junction map
            auto& junction = junctions[key];
            junction.add_sample(sample_name, 1);
        }

        // Update reference position for operations that consume reference
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
            op == BAM_CEQUAL || op == BAM_CDIFF) {
            ref_pos += len;
        }
    }
}

bool ParallelBamExtractor::has_bam_index(const std::string& bam_path) {
    // Check for .bai or .bam.bai
    std::string bai_path1 = bam_path + ".bai";
    std::string bai_path2 = bam_path.substr(0, bam_path.length() - 4) + ".bai";

    return fs::exists(bai_path1) || fs::exists(bai_path2);
}

JunctionMap ParallelBamExtractor::extract_junctions_from_bam(
    const std::string& bam_path,
    const std::string& sample_name,
    int htslib_threads) {

    JunctionMap junctions;

    // Open BAM file
    samFile* in = sam_open(bam_path.c_str(), "r");
    if (!in) {
        throw std::runtime_error("Cannot open BAM file: " + bam_path);
    }

    // Enable multi-threaded BGZF decompression
    if (htslib_threads > 1) {
        hts_set_threads(in, htslib_threads);
    }

    // Read header
    sam_hdr_t* header = sam_hdr_read(in);
    if (!header) {
        sam_close(in);
        throw std::runtime_error("Cannot read BAM header: " + bam_path);
    }

    // Allocate alignment record
    bam1_t* aln = bam_init1();

    // Read all alignments
    while (sam_read1(in, header, aln) >= 0) {
        extract_junctions_from_alignment(aln, header, sample_name, junctions);
    }

    // Cleanup
    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(in);

    return junctions;
}

JunctionMap ParallelBamExtractor::extract_junctions_from_region(
    const std::string& bam_path,
    const std::string& sample_name,
    const std::string& region) {

    JunctionMap junctions;

    // Open BAM file
    samFile* in = sam_open(bam_path.c_str(), "r");
    if (!in) {
        throw std::runtime_error("Cannot open BAM file: " + bam_path);
    }

    // Read header
    sam_hdr_t* header = sam_hdr_read(in);
    if (!header) {
        sam_close(in);
        throw std::runtime_error("Cannot read BAM header: " + bam_path);
    }

    // Load index
    hts_idx_t* idx = sam_index_load(in, bam_path.c_str());
    if (!idx) {
        sam_hdr_destroy(header);
        sam_close(in);
        throw std::runtime_error("Cannot load BAM index for: " + bam_path);
    }

    // Create iterator for region
    hts_itr_t* iter = sam_itr_querys(idx, header, region.c_str());
    if (!iter) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(header);
        sam_close(in);
        // Region might not exist in this BAM, return empty
        return junctions;
    }

    // Allocate alignment record
    bam1_t* aln = bam_init1();

    // Read alignments in region
    while (sam_itr_next(in, iter, aln) >= 0) {
        extract_junctions_from_alignment(aln, header, sample_name, junctions);
    }

    // Cleanup
    bam_destroy1(aln);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(in);

    return junctions;
}

JunctionMap ParallelBamExtractor::extract_junctions_from_bam_tiled(
    const std::string& bam_path,
    const std::string& sample_name,
    int num_threads,
    int64_t tile_size,
    bool verbose) {

    // Check for index
    if (!has_bam_index(bam_path)) {
        if (verbose) {
            std::cerr << "  No BAM index found, falling back to sequential processing\n";
        }
        return extract_junctions_from_bam(bam_path, sample_name, num_threads);
    }

    // Open BAM to get header info
    samFile* in = sam_open(bam_path.c_str(), "r");
    if (!in) {
        throw std::runtime_error("Cannot open BAM file: " + bam_path);
    }

    sam_hdr_t* header = sam_hdr_read(in);
    if (!header) {
        sam_close(in);
        throw std::runtime_error("Cannot read BAM header: " + bam_path);
    }

    // Build list of regions (tiles)
    std::vector<std::string> regions;
    int n_targets = sam_hdr_nref(header);

    for (int tid = 0; tid < n_targets; ++tid) {
        const char* chrom = sam_hdr_tid2name(header, tid);
        int64_t chrom_len = sam_hdr_tid2len(header, tid);

        // Create tiles for this chromosome
        for (int64_t start = 0; start < chrom_len; start += tile_size) {
            int64_t end = std::min(start + tile_size, chrom_len);
            std::ostringstream region;
            region << chrom << ":" << (start + 1) << "-" << end;  // 1-based for samtools
            regions.push_back(region.str());
        }
    }

    sam_hdr_destroy(header);
    sam_close(in);

    if (verbose) {
        std::cerr << "  Processing " << regions.size() << " tiles with "
                  << num_threads << " threads\n";
    }

    // Process regions in parallel
    ThreadPool pool(num_threads);
    std::vector<std::future<JunctionMap>> futures;

    for (const auto& region : regions) {
        futures.push_back(pool.submit([&bam_path, &sample_name, region]() {
            return extract_junctions_from_region(bam_path, sample_name, region);
        }));
    }

    // Collect results
    std::vector<JunctionMap> results;
    results.reserve(futures.size());
    for (auto& future : futures) {
        results.push_back(future.get());
    }

    // Merge all region results
    return merge_junction_maps(results);
}

JunctionMap ParallelBamExtractor::merge_junction_maps(std::vector<JunctionMap>& maps) {
    if (maps.empty()) {
        return JunctionMap();
    }

    // Use the largest map as the base to minimize copying
    size_t largest_idx = 0;
    size_t largest_size = maps[0].size();
    for (size_t i = 1; i < maps.size(); ++i) {
        if (maps[i].size() > largest_size) {
            largest_size = maps[i].size();
            largest_idx = i;
        }
    }

    JunctionMap merged = std::move(maps[largest_idx]);

    // Merge other maps into the base
    for (size_t i = 0; i < maps.size(); ++i) {
        if (i == largest_idx) continue;

        for (auto& [key, data] : maps[i]) {
            auto it = merged.find(key);
            if (it == merged.end()) {
                merged[key] = std::move(data);
            } else {
                // Combine scores from this sample
                for (const auto& sample : data.samples) {
                    it->second.add_sample(sample.sample_name, sample.score);
                }
            }
        }
    }

    return merged;
}

JunctionMap ParallelBamExtractor::extract_junctions_parallel(
    const std::vector<std::string>& bam_paths,
    int num_threads,
    bool verbose) {

    if (bam_paths.empty()) {
        return JunctionMap();
    }

    // Get sample names from file paths
    auto get_sample_name = [](const std::string& path) {
        return fs::path(path).stem().string();
    };

    size_t num_files = bam_paths.size();

    // Decide parallelism strategy
    // If we have more threads than files, use tiled processing within files
    bool use_tiled = (num_threads > static_cast<int>(num_files)) &&
                     std::all_of(bam_paths.begin(), bam_paths.end(), has_bam_index);

    if (verbose) {
        std::cerr << "Extracting junctions from " << num_files << " BAM files using "
                  << num_threads << " threads";
        if (use_tiled) {
            std::cerr << " (tiled intra-file parallelism)";
        } else {
            std::cerr << " (inter-file parallelism + htslib threading)";
        }
        std::cerr << "...\n";
    }

    // Progress tracking
    std::atomic<size_t> completed_files{0};
    std::mutex progress_mutex;

    std::vector<JunctionMap> results;
    results.resize(num_files);

    if (use_tiled) {
        // Tiled processing: distribute threads across files
        // Each file gets approximately equal thread allocation
        int threads_per_file = std::max(1, num_threads / static_cast<int>(num_files));
        int extra_threads = num_threads - (threads_per_file * static_cast<int>(num_files));

        std::vector<std::future<JunctionMap>> futures;

        for (size_t i = 0; i < num_files; ++i) {
            const std::string& bam_path = bam_paths[i];
            std::string sample_name = get_sample_name(bam_path);

            // Give extra threads to first few files
            int file_threads = threads_per_file + (static_cast<int>(i) < extra_threads ? 1 : 0);

            futures.push_back(std::async(std::launch::async,
                [&, i, bam_path, sample_name, file_threads]() {
                    if (verbose) {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cerr << "  [" << (i + 1) << "/" << num_files
                                  << "] Starting: " << sample_name
                                  << " (" << file_threads << " threads)\n";
                    }

                    JunctionMap result = extract_junctions_from_bam_tiled(
                        bam_path, sample_name, file_threads, 10000000, false);

                    size_t done = ++completed_files;
                    if (verbose) {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cerr << "  [" << done << "/" << num_files
                                  << "] Completed: " << sample_name
                                  << " (" << result.size() << " junctions)\n";
                    }

                    return result;
                }));
        }

        for (size_t i = 0; i < futures.size(); ++i) {
            results[i] = futures[i].get();
        }
    } else {
        // Inter-file parallelism with htslib multi-threaded decompression
        // Calculate htslib threads per file
        int htslib_threads = std::max(1, num_threads / static_cast<int>(num_files));

        auto process_bam = [&](size_t idx) -> JunctionMap {
            const std::string& bam_path = bam_paths[idx];
            std::string sample_name = get_sample_name(bam_path);

            if (verbose) {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cerr << "  [" << (idx + 1) << "/" << num_files
                          << "] Starting: " << sample_name << "\n";
            }

            JunctionMap result = extract_junctions_from_bam(bam_path, sample_name, htslib_threads);

            size_t done = ++completed_files;
            if (verbose) {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cerr << "  [" << done << "/" << num_files
                          << "] Completed: " << sample_name
                          << " (" << result.size() << " junctions)\n";
            }

            return result;
        };

        // Launch extraction tasks
        std::vector<std::future<JunctionMap>> futures;

        for (size_t i = 0; i < num_files; ++i) {
            futures.push_back(std::async(
                num_threads > 1 ? std::launch::async : std::launch::deferred,
                process_bam, i));
        }

        // Wait for all tasks and collect results
        for (size_t i = 0; i < futures.size(); ++i) {
            results[i] = futures[i].get();
        }
    }

    if (verbose) {
        std::cerr << "Merging junction maps...\n";
    }

    // Merge all results
    JunctionMap merged = merge_junction_maps(results);

    if (verbose) {
        std::cerr << "Total unique junctions: " << merged.size() << "\n";
    }

    return merged;
}

} // namespace eastr
