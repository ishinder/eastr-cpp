#include "eastr/types.hpp"
#include "eastr/junction.hpp"
#include "eastr/bed_parser.hpp"
#include "eastr/gtf_parser.hpp"
#include "eastr/fasta_index.hpp"
#include "eastr/junction_extractor.hpp"
#include "eastr/spurious_detector.hpp"
#include "eastr/bam_filter.hpp"
#include "eastr/output_writer.hpp"
#include "eastr/multi_bam_reader.hpp"
#include "eastr/junction_accumulator.hpp"
#include "eastr/streaming_spurious_detector.hpp"
#include "eastr/parallel_bam_extractor.hpp"

#include <CLI/CLI.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <tuple>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

// Check if bowtie2 index exists
bool bowtie2_index_exists(const std::string& index_prefix) {
    // Check for .1.bt2 or .1.bt2l (large index)
    return fs::exists(index_prefix + ".1.bt2") || fs::exists(index_prefix + ".1.bt2l");
}

// Get base name from FASTA path (removes .fa/.fasta/.gz extensions)
std::string get_fasta_basename(const std::string& fasta_path) {
    std::string base_name = fs::path(fasta_path).stem().string();

    // Remove .fa or .fasta from stem if present (for .fa.gz files)
    if (base_name.size() >= 3 && base_name.substr(base_name.size() - 3) == ".fa") {
        base_name = base_name.substr(0, base_name.size() - 3);
    } else if (base_name.size() >= 6 && base_name.substr(base_name.size() - 6) == ".fasta") {
        base_name = base_name.substr(0, base_name.size() - 6);
    }
    return base_name;
}

// Build bowtie2 index from FASTA with lock file to prevent concurrent builds
// Returns: tuple of (index_prefix, index_dir, created_by_eastr)
// - index_prefix: path to use with bowtie2
// - index_dir: directory containing the index
// - created_by_eastr: true if this process built the index
std::tuple<std::string, std::string, bool> build_bowtie2_index(
    const std::string& fasta_path,
    const std::string& output_dir,
    bool verbose) {

    std::string base_name = get_fasta_basename(fasta_path);
    std::string index_dir_name = base_name + "_bt2_idx";

    // Check if index already exists alongside the FASTA
    std::string fasta_dir_prefix = (fs::path(fasta_path).parent_path() / base_name).string();
    if (bowtie2_index_exists(fasta_dir_prefix)) {
        if (verbose) {
            std::cerr << "Using existing bowtie2 index: " << fasta_dir_prefix << "\n";
        }
        return {fasta_dir_prefix, "", false};
    }

    // Determine where to put the index
    fs::path index_dir;
    if (!output_dir.empty()) {
        index_dir = fs::path(output_dir) / index_dir_name;
    } else {
        // Use directory alongside FASTA if no output specified
        index_dir = fs::path(fasta_path).parent_path() / index_dir_name;
    }

    std::string index_prefix = (index_dir / base_name).string();

    // Check if index already exists in output directory
    if (bowtie2_index_exists(index_prefix)) {
        if (verbose) {
            std::cerr << "Using existing bowtie2 index: " << index_prefix << "\n";
        }
        return {index_prefix, "", false};
    }

    // Create directory for index and lock file
    fs::create_directories(index_dir);

    // Use a lock file to prevent concurrent builds
    std::string lock_file = (index_dir / ".building.lock").string();

    // Try to acquire lock (O_CREAT | O_EXCL fails if file exists)
    int lock_fd = open(lock_file.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);

    if (lock_fd >= 0) {
        // We acquired the lock - we're responsible for building
        close(lock_fd);

        // Double-check index wasn't created between our check and lock acquisition
        if (bowtie2_index_exists(index_prefix)) {
            fs::remove(lock_file);
            if (verbose) {
                std::cerr << "Using existing bowtie2 index: " << index_prefix << "\n";
            }
            return {index_prefix, "", false};
        }

        if (verbose) {
            std::cerr << "Building bowtie2 index from " << fasta_path << "...\n";
            std::cerr << "  Index directory: " << index_dir.string() << "\n";
        }

        // Build the index
        std::string cmd = "bowtie2-build";
        cmd += " \"" + fasta_path + "\"";
        cmd += " \"" + index_prefix + "\"";

        if (!verbose) {
            cmd += " >/dev/null 2>&1";
        }

        int ret = std::system(cmd.c_str());

        // Remove lock file (whether success or failure)
        fs::remove(lock_file);

        if (ret != 0) {
            fs::remove_all(index_dir);
            throw std::runtime_error("Failed to build bowtie2 index. Make sure bowtie2-build is installed and in PATH.");
        }

        if (verbose) {
            std::cerr << "Bowtie2 index built: " << index_prefix << "\n";
        }

        return {index_prefix, index_dir.string(), true};

    } else {
        // Lock file exists - another process is building the index
        if (verbose) {
            std::cerr << "Waiting for another process to finish building bowtie2 index...\n";
        }

        // Wait for lock file to disappear (building to complete)
        int wait_seconds = 0;
        const int max_wait = 3600;  // Max 1 hour wait
        while (fs::exists(lock_file) && wait_seconds < max_wait) {
            std::this_thread::sleep_for(std::chrono::seconds(5));
            wait_seconds += 5;
            if (verbose && wait_seconds % 60 == 0) {
                std::cerr << "  Still waiting... (" << wait_seconds << "s)\n";
            }
        }

        if (wait_seconds >= max_wait) {
            throw std::runtime_error("Timeout waiting for bowtie2 index to be built by another process. "
                "Lock file: " + lock_file);
        }

        // Verify index was actually built
        if (!bowtie2_index_exists(index_prefix)) {
            throw std::runtime_error("Another process was building bowtie2 index but it failed. "
                "Please check for errors or remove directory: " + index_dir.string());
        }

        if (verbose) {
            std::cerr << "Using bowtie2 index built by another process: " << index_prefix << "\n";
        }

        return {index_prefix, index_dir.string(), false};
    }
}

// Check if path is a file (has extension) or directory
bool is_file_path(const std::string& path) {
    size_t last_dot = path.find_last_of('.');
    size_t last_slash = path.find_last_of("/\\");
    return last_dot != std::string::npos &&
           (last_slash == std::string::npos || last_dot > last_slash);
}

// Expand file list (if .txt file, read paths from it)
std::vector<std::string> expand_file_list(const std::string& path) {
    std::vector<std::string> files;

    // Check extension
    size_t dot_pos = path.find_last_of('.');
    std::string ext = (dot_pos != std::string::npos) ? path.substr(dot_pos) : "";

    if (ext == ".bam" || ext == ".cram" || ext == ".sam" ||
        ext == ".bed" || ext == ".gtf" || ext == ".gff") {
        // Single file
        files.push_back(path);
    } else {
        // Assume it's a list file
        std::ifstream list_file(path);
        if (!list_file.is_open()) {
            throw std::runtime_error("Cannot open file list: " + path);
        }

        std::string line;
        while (std::getline(list_file, line)) {
            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);

            if (!line.empty()) {
                if (!fs::exists(line)) {
                    throw std::runtime_error("File not found: " + line);
                }
                files.push_back(line);
            }
        }
    }

    return files;
}

int main(int argc, char* argv[]) {
    CLI::App app{"eastr: Emending alignments of spuriously spliced transcript reads"};

    eastr::Config config;

    // Input options (mutually exclusive)
    auto* gtf_opt = app.add_option("--gtf", config.gtf_path,
        "Input GTF file containing transcript annotations");
    auto* bed_opt = app.add_option("--bed", config.bed_path,
        "Input BED file with intron coordinates");
    auto* bam_opt = app.add_option("--bam", config.bam_path,
        "Input BAM file or TXT file containing list of BAM files");

    gtf_opt->excludes(bed_opt)->excludes(bam_opt);
    bed_opt->excludes(gtf_opt)->excludes(bam_opt);
    bam_opt->excludes(gtf_opt)->excludes(bed_opt);

    // Required options
    app.add_option("-r,--reference", config.reference_fasta,
        "Reference FASTA genome (uncompressed or bgzip-compressed)")->required();
    app.add_option("-i,--bowtie2_index", config.bowtie2_index,
        "Path to Bowtie2 index prefix (optional, will be built from reference if not provided)");

    // Algorithm parameters
    app.add_option("--bt2_k", config.params.bt2_k,
        "Min distinct alignments for spurious classification (default: 10)");
    app.add_option("-o", config.params.overhang,
        "Length of overhang on either side of splice junction (default: 50)");
    app.add_option("-a", config.params.anchor,
        "Minimum required anchor length in each exon (default: 7)");
    app.add_option("--min_duplicate_exon_length", config.params.min_duplicate_exon_length,
        "Minimum length for duplicated exon detection (default: 27)");
    app.add_option("--min_junc_score", config.params.min_junc_score,
        "Minimum supporting reads per junction (default: 1)");

    // Minimap2/alignment parameters
    app.add_option("-A", config.scoring.match_score,
        "Matching score (default: 3)");
    app.add_option("-B", config.scoring.mismatch_penalty,
        "Mismatching penalty (default: 4)");
    app.add_option("-k", config.params.kmer_size,
        "K-mer length for alignment (default: 3)");
    app.add_option("-w", config.params.window_size,
        "Minimizer window size (default: 2)");
    app.add_option("-m", config.params.min_chain_score,
        "Minimum chain score (default: 25)");

    // Output options
    app.add_option("--out_original_junctions", config.out_original_junctions,
        "Write original junctions to file/directory");
    app.add_option("--out_removed_junctions", config.out_removed_junctions,
        "Write spurious junctions to file (default: stdout)");
    app.add_option("--out_kept_junctions", config.out_kept_junctions,
        "Write non-spurious junctions to file (optional)");
    app.add_option("--out_filtered_bam", config.out_filtered_bam,
        "Write filtered BAMs to file/directory");
    app.add_option("--filtered_bam_suffix", config.filtered_bam_suffix,
        "Suffix for filtered BAM files (default: _EASTR_filtered)");

    // Other options
    app.add_option("-p", config.num_threads,
        "Number of parallel processes (default: 1)");
    app.add_option("--batch-size", config.batch_size,
        "Number of junctions to batch for spurious testing (default: 1000)");
    app.add_option("--prefetch-buffer", config.prefetch_buffer_size,
        "Prefetch buffer size for BAM reading (default: 10000)");
    app.add_flag("--verbose", config.verbose,
        "Display additional information during processing");
    app.add_flag("--removed_alignments_bam", config.write_removed_alignments,
        "Write removed alignments to a BAM file");
    app.add_option("--trusted_bed", config.trusted_bed,
        "BED file with trusted junctions (will not be removed)");

    CLI11_PARSE(app, argc, argv);

    // Validate input
    if (config.gtf_path.empty() && config.bed_path.empty() && config.bam_path.empty()) {
        std::cerr << "Error: One of --gtf, --bed, or --bam is required\n";
        return 1;
    }

    // Track bt2 index for potential cleanup (only if EASTR created it)
    std::string bt2_index_dir;
    bool bt2_created_by_eastr = false;

    try {
        // Index reference if needed
        eastr::FastaIndex::ensure_indexed(config.reference_fasta);

        // Determine output directory for index
        std::string output_dir;
        if (!config.out_filtered_bam.empty()) {
            output_dir = fs::path(config.out_filtered_bam).parent_path().string();
            if (output_dir.empty()) output_dir = config.out_filtered_bam;
        } else if (!config.out_removed_junctions.empty() && config.out_removed_junctions != "stdout") {
            output_dir = fs::path(config.out_removed_junctions).parent_path().string();
        }

        // Build or locate bowtie2 index
        std::string bt2_index = config.bowtie2_index;
        if (bt2_index.empty()) {
            auto [index_path, index_dir, created] = build_bowtie2_index(config.reference_fasta, output_dir, config.verbose);
            bt2_index = index_path;
            bt2_index_dir = index_dir;
            bt2_created_by_eastr = created;
        } else if (!bowtie2_index_exists(bt2_index)) {
            throw std::runtime_error("Bowtie2 index not found: " + bt2_index +
                "\nEither provide a valid index path or omit -i to build automatically.");
        }

        // Determine input type and extract junctions
        eastr::JunctionMap junctions;
        eastr::JunctionMap spurious;  // Declare early for streaming path
        bool is_bam_input = false;
        std::vector<std::string> bam_files;
        std::vector<std::string> sample_names;  // For consistent output ordering

        if (!config.bam_path.empty()) {
            is_bam_input = true;
            bam_files = expand_file_list(config.bam_path);

            // Get sample names (for consistent output ordering)
            for (size_t i = 0; i < bam_files.size(); ++i) {
                fs::path p(bam_files[i]);
                sample_names.push_back(p.stem().string());
            }

            if (config.verbose) {
                std::cerr << "Extracting junctions from BAM files in parallel...\n";
                std::cerr << "  Files: " << bam_files.size() << "\n";
                std::cerr << "  Threads: " << config.num_threads << "\n";
            }

            // Extract junctions from all BAM files in parallel
            junctions = eastr::ParallelBamExtractor::extract_junctions_parallel(
                bam_files, config.num_threads, config.verbose);

            // Write original junctions if requested
            if (!config.out_original_junctions.empty()) {
                std::vector<std::string> orig_junction_paths =
                    eastr::OutputWriter::generate_junction_output_paths(
                        bam_files, config.out_original_junctions, "_original_junctions");

                eastr::OutputWriter::write_spurious_bed(
                    junctions, config.scoring, orig_junction_paths[0],
                    eastr::OutputWriter::InputType::BAM, sample_names);
            }

            if (config.verbose) {
                std::cerr << "Total unique junctions: " << junctions.size() << "\n";
                std::cerr << "Detecting spurious junctions...\n";
            }

            // Detect spurious junctions using batch approach
            eastr::SpuriousDetector batch_detector(
                config.scoring,
                config.params,
                bt2_index,
                config.reference_fasta,
                config.num_threads);

            std::optional<std::string> trusted_bed;
            if (!config.trusted_bed.empty()) {
                trusted_bed = config.trusted_bed;
            }

            spurious = batch_detector.detect_spurious(
                junctions, is_bam_input, trusted_bed);

            if (config.verbose) {
                std::cerr << "Spurious junctions: " << spurious.size() << "\n";
            }

        } else if (!config.bed_path.empty()) {
            auto bed_files = expand_file_list(config.bed_path);

            if (config.verbose) {
                std::cerr << "Parsing BED files...\n";
            }

            junctions = eastr::BedParser::parse_multi(bed_files, config.num_threads);

        } else if (!config.gtf_path.empty()) {
            if (config.verbose) {
                std::cerr << "Parsing GTF file...\n";
            }

            junctions = eastr::GtfParser::parse(config.gtf_path);
        }

        // For non-BAM inputs, use batch spurious detection
        // (BAM inputs already populated spurious via streaming)
        if (!is_bam_input) {
            if (config.verbose) {
                std::cerr << "Detecting spurious junctions...\n";
            }

            // Detect spurious junctions using batch approach
            eastr::SpuriousDetector batch_detector(
                config.scoring,
                config.params,
                bt2_index,
                config.reference_fasta,
                config.num_threads);

            std::optional<std::string> trusted_bed;
            if (!config.trusted_bed.empty()) {
                trusted_bed = config.trusted_bed;
            }

            spurious = batch_detector.detect_spurious(
                junctions, is_bam_input, trusted_bed);

            if (config.verbose) {
                std::cerr << "Found " << spurious.size() << " spurious junctions\n";
            }
        }

        // Determine input type for output formatting
        eastr::OutputWriter::InputType output_type;
        if (!config.gtf_path.empty()) {
            output_type = eastr::OutputWriter::InputType::GTF;
        } else if (!config.bed_path.empty()) {
            output_type = eastr::OutputWriter::InputType::BED;
        } else {
            output_type = eastr::OutputWriter::InputType::BAM;
        }

        // Output results
        if (is_bam_input && !config.out_filtered_bam.empty()) {
            // Write spurious junctions to a single BED file (used for all BAM filtering)
            std::string spurious_bed_path = config.out_removed_junctions;

            // If output is "stdout", we need a temp file for BAM filtering
            bool using_temp_bed = (spurious_bed_path == "stdout" || spurious_bed_path.empty());
            if (using_temp_bed) {
                // Create temp file in the filtered BAM output directory
                fs::path temp_dir = fs::path(config.out_filtered_bam).parent_path();
                if (temp_dir.empty()) temp_dir = config.out_filtered_bam;
                fs::create_directories(temp_dir);
                spurious_bed_path = (temp_dir / ".eastr_spurious_temp.bed").string();
            } else {
                // Create output directory if needed
                fs::path dir = fs::path(spurious_bed_path).parent_path();
                if (!dir.empty()) {
                    fs::create_directories(dir);
                }
            }

            // Write the combined spurious junctions BED file
            eastr::OutputWriter::write_spurious_bed(
                spurious, config.scoring, spurious_bed_path, output_type, sample_names);

            // Generate filtered BAM output paths
            auto filtered_bam_paths = eastr::OutputWriter::generate_bam_output_paths(
                bam_files, config.out_filtered_bam, config.filtered_bam_suffix);

            // Create output directories if needed
            for (const auto& path : filtered_bam_paths) {
                fs::path dir = fs::path(path).parent_path();
                if (!dir.empty()) {
                    fs::create_directories(dir);
                }
            }

            if (config.verbose) {
                std::cerr << "Filtering BAM files...\n";
                int io_threads_per_file = std::max(1, config.num_threads / static_cast<int>(bam_files.size()));
                if (io_threads_per_file > 1) {
                    std::cerr << "  Using " << io_threads_per_file << " I/O threads per file for BGZF compression\n";
                }
            }

            eastr::BamFilter::Options filter_opts;
            filter_opts.remove_mate = true;
            filter_opts.verbose = config.verbose;
            filter_opts.write_removed = config.write_removed_alignments;

            // Filter all BAMs using the single spurious junctions BED file
            eastr::filter_multi_bam(
                bam_files, spurious_bed_path, filtered_bam_paths,
                config.num_threads, filter_opts);

            // Clean up temp file if we created one
            if (using_temp_bed) {
                fs::remove(spurious_bed_path);
            }

        } else {
            // Just output the spurious junctions
            eastr::OutputWriter::write_spurious_bed(
                spurious, config.scoring, config.out_removed_junctions, output_type, sample_names);
        }

        // Write kept (non-spurious) junctions if requested
        if (!config.out_kept_junctions.empty()) {
            // Build map of kept junctions (all junctions minus spurious)
            eastr::JunctionMap kept;
            for (const auto& [key, data] : junctions) {
                if (spurious.find(key) == spurious.end()) {
                    kept[key] = data;
                }
            }

            if (config.verbose) {
                std::cerr << "Writing " << kept.size() << " non-spurious junctions to "
                          << config.out_kept_junctions << "\n";
            }

            // Create output directory if needed
            fs::path dir = fs::path(config.out_kept_junctions).parent_path();
            if (!dir.empty()) {
                fs::create_directories(dir);
            }

            eastr::OutputWriter::write_spurious_bed(
                kept, config.scoring, config.out_kept_junctions, output_type, sample_names);
        }

        // Note: Auto-built bowtie2 index is preserved for reuse by concurrent/future runs
        // Users can manually delete it when all jobs are complete
        if (bt2_created_by_eastr && config.verbose) {
            std::cerr << "Bowtie2 index saved at: " << bt2_index_dir << "\n";
        }

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
