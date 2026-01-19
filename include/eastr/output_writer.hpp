#pragma once

#include "junction.hpp"
#include "types.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace eastr {

class OutputWriter {
public:
    enum class InputType { BAM, GTF, BED };

    // Write spurious junctions to BED format
    // sample_order: optional list of sample names to define output order
    static void write_spurious_bed(
        const JunctionMap& spurious,
        const ScoringParams& scoring,
        const std::string& output_path,  // "stdout" for console output
        InputType input_type,
        const std::vector<std::string>& sample_order = {});

    // Write junctions grouped by sample
    // Returns map of sample_name -> bed_file_path
    static std::unordered_map<std::string, std::string>
    write_spurious_by_sample(
        const JunctionMap& spurious,
        const std::vector<std::string>& bam_paths,
        const std::vector<std::string>& output_paths,
        const ScoringParams& scoring);

    // Generate output file paths based on input
    static std::vector<std::string> generate_junction_output_paths(
        const std::vector<std::string>& input_files,
        const std::string& output_dir_or_file,
        const std::string& suffix);

    static std::vector<std::string> generate_bam_output_paths(
        const std::vector<std::string>& input_bams,
        const std::string& output_dir_or_file,
        const std::string& suffix);
};

} // namespace eastr
