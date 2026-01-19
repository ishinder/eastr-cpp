#include "eastr/output_writer.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace eastr {

void OutputWriter::write_spurious_bed(
    const JunctionMap& spurious,
    const ScoringParams& scoring,
    const std::string& output_path,
    InputType input_type,
    const std::vector<std::string>& sample_order) {

    // Build sample name -> index map for sorting
    std::unordered_map<std::string, size_t> sample_index;
    for (size_t i = 0; i < sample_order.size(); ++i) {
        sample_index[sample_order[i]] = i;
    }

    // Helper to sort samples by their order in sample_order
    auto sort_samples = [&sample_index](std::vector<SampleInfo> samples) {
        if (sample_index.empty()) {
            return samples;  // No ordering specified, keep as-is
        }
        std::sort(samples.begin(), samples.end(),
            [&sample_index](const SampleInfo& a, const SampleInfo& b) {
                auto it_a = sample_index.find(a.sample_name);
                auto it_b = sample_index.find(b.sample_name);
                size_t idx_a = (it_a != sample_index.end()) ? it_a->second : SIZE_MAX;
                size_t idx_b = (it_b != sample_index.end()) ? it_b->second : SIZE_MAX;
                return idx_a < idx_b;
            });
        return samples;
    };

    // Sort junctions for consistent output
    std::vector<std::pair<JunctionKey, const JunctionData*>> sorted;
    for (const auto& [key, data] : spurious) {
        sorted.emplace_back(key, &data);
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    // Determine output stream
    std::ostream* out_stream = nullptr;
    std::ofstream file_stream;

    if (output_path == "stdout" || output_path.empty()) {
        out_stream = &std::cout;
    } else {
        file_stream.open(output_path);
        if (!file_stream.is_open()) {
            throw std::runtime_error("Cannot open output file: " + output_path);
        }
        out_stream = &file_stream;
    }

    // Determine number of digits for junction numbering
    int num_digits = std::to_string(sorted.size()).size();

    int junction_num = 1;
    for (const auto& [key, data] : sorted) {
        std::ostringstream name;
        name << "JUNC" << std::setw(num_digits) << std::setfill('0') << junction_num++;

        // Calculate alignment score if hit exists
        int align_score = data->hit ? data->hit->score : 0;

        *out_stream << key.chrom << "\t"
                    << key.start << "\t"
                    << key.end << "\t"
                    << name.str() << "\t";

        if (input_type == InputType::GTF) {
            // GTF format: score is '.', include gene/transcript info
            *out_stream << "." << "\t"
                        << to_char(key.strand) << "\t"
                        << align_score << "\t";

            // Gene and transcript IDs
            *out_stream << data->gene_id;
            for (const auto& tid : data->transcript_ids) {
                *out_stream << ";" << tid;
            }
        } else if (input_type == InputType::BED) {
            // BED format: include sample info
            *out_stream << data->total_score << "\t"
                        << to_char(key.strand) << "\t";

            // Sample info (sorted by input order)
            auto sorted_samples = sort_samples(data->samples);
            bool first = true;
            for (const auto& sample : sorted_samples) {
                if (!first) *out_stream << ";";
                *out_stream << sample.sample_name << "," << sample.score;
                first = false;
            }
        } else {
            // BAM format
            *out_stream << data->total_score << "\t"
                        << to_char(key.strand) << "\t"
                        << align_score << "\t";

            // Sample info (sorted by input order)
            auto sorted_samples = sort_samples(data->samples);
            bool first = true;
            for (const auto& sample : sorted_samples) {
                if (!first) *out_stream << ";";
                *out_stream << sample.score << "," << sample.sample_name;
                first = false;
            }
        }

        *out_stream << "\n";
    }
}

std::unordered_map<std::string, std::string>
OutputWriter::write_spurious_by_sample(
    const JunctionMap& spurious,
    const std::vector<std::string>& bam_paths,
    const std::vector<std::string>& output_paths,
    const ScoringParams& scoring) {

    // Get sample names from BAM paths
    auto get_sample_name = [](const std::string& path) {
        size_t last_slash = path.find_last_of("/\\");
        size_t last_dot = path.find_last_of('.');
        if (last_slash == std::string::npos) last_slash = 0;
        else last_slash++;
        if (last_dot == std::string::npos || last_dot < last_slash) {
            return path.substr(last_slash);
        }
        return path.substr(last_slash, last_dot - last_slash);
    };

    // Create sample name -> output path mapping
    std::unordered_map<std::string, std::string> sample_to_bed;
    std::unordered_map<std::string, std::ofstream> sample_files;

    for (size_t i = 0; i < bam_paths.size(); ++i) {
        std::string sample = get_sample_name(bam_paths[i]);
        std::string output = (i < output_paths.size()) ? output_paths[i] : "";

        if (!output.empty()) {
            sample_to_bed[sample] = output;
            sample_files[sample].open(output);
        }
    }

    // Sort junctions
    std::vector<std::pair<JunctionKey, const JunctionData*>> sorted;
    for (const auto& [key, data] : spurious) {
        sorted.emplace_back(key, &data);
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    int num_digits = std::to_string(sorted.size()).size();

    // Write junctions to per-sample files
    int junction_num = 1;
    for (const auto& [key, data] : sorted) {
        std::ostringstream name;
        name << "JUNC" << std::setw(num_digits) << std::setfill('0') << junction_num++;

        int align_score = data->hit ? data->hit->score : 0;

        for (const auto& sample_info : data->samples) {
            auto file_it = sample_files.find(sample_info.sample_name);
            if (file_it == sample_files.end() || !file_it->second.is_open()) {
                continue;
            }

            file_it->second << key.chrom << "\t"
                            << key.start << "\t"
                            << key.end << "\t"
                            << name.str() << "\t"
                            << sample_info.score << "\t"
                            << to_char(key.strand) << "\t"
                            << align_score << "\n";
        }
    }

    // Close all files
    for (auto& [sample, file] : sample_files) {
        file.close();
    }

    return sample_to_bed;
}

std::vector<std::string> OutputWriter::generate_junction_output_paths(
    const std::vector<std::string>& input_files,
    const std::string& output_dir_or_file,
    const std::string& suffix) {

    std::vector<std::string> output_paths;

    auto get_basename = [](const std::string& path) {
        size_t last_slash = path.find_last_of("/\\");
        size_t last_dot = path.find_last_of('.');
        if (last_slash == std::string::npos) last_slash = 0;
        else last_slash++;
        if (last_dot == std::string::npos || last_dot < last_slash) {
            return path.substr(last_slash);
        }
        return path.substr(last_slash, last_dot - last_slash);
    };

    auto get_parent_dir = [](const std::string& path) {
        size_t last_slash = path.find_last_of("/\\");
        if (last_slash == std::string::npos) {
            return std::string(".");
        }
        return path.substr(0, last_slash);
    };

    // Check if output is a file or directory
    bool is_file = output_dir_or_file.find('.') != std::string::npos &&
                   output_dir_or_file.find_last_of('.') > output_dir_or_file.find_last_of("/\\");

    if (input_files.size() == 1 && is_file) {
        output_paths.push_back(output_dir_or_file);
    } else {
        // Determine the output directory
        // If output looks like a file but we have multiple inputs, use parent directory
        std::string out_dir = is_file ? get_parent_dir(output_dir_or_file) : output_dir_or_file;

        for (const auto& input : input_files) {
            std::string basename = get_basename(input);
            output_paths.push_back(out_dir + "/" + basename + suffix + ".bed");
        }
    }

    return output_paths;
}

std::vector<std::string> OutputWriter::generate_bam_output_paths(
    const std::vector<std::string>& input_bams,
    const std::string& output_dir_or_file,
    const std::string& suffix) {

    std::vector<std::string> output_paths;

    auto get_basename = [](const std::string& path) {
        size_t last_slash = path.find_last_of("/\\");
        size_t last_dot = path.find_last_of('.');
        if (last_slash == std::string::npos) last_slash = 0;
        else last_slash++;
        if (last_dot == std::string::npos || last_dot < last_slash) {
            return path.substr(last_slash);
        }
        return path.substr(last_slash, last_dot - last_slash);
    };

    bool is_file = output_dir_or_file.find('.') != std::string::npos &&
                   output_dir_or_file.find_last_of('.') > output_dir_or_file.find_last_of("/\\");

    if (input_bams.size() == 1 && is_file) {
        output_paths.push_back(output_dir_or_file);
    } else {
        for (const auto& input : input_bams) {
            std::string basename = get_basename(input);
            output_paths.push_back(output_dir_or_file + "/" + basename + suffix + ".bam");
        }
    }

    return output_paths;
}

} // namespace eastr
