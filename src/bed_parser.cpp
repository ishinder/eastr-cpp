#include "eastr/bed_parser.hpp"

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <future>
#include <sstream>
#include <stdexcept>

namespace eastr {

// Helper to check if file is gzipped
static bool is_gzipped(const std::string& path) {
    return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
}

// Process a single BED line
static bool process_bed_line(const std::string& line, const std::string& bed_path,
                             std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash>& junctions) {
    // Skip empty lines
    if (line.empty()) return true;

    // Skip comment lines and track lines
    if (line[0] == '#' || line.substr(0, 5) == "track") return true;

    std::istringstream iss(line);
    std::string chrom, name, strand_str;
    int64_t start, end;
    int score;

    // Parse BED6 format: chrom start end name score strand
    if (!(iss >> chrom >> start >> end >> name >> score >> strand_str)) {
        throw std::runtime_error("Invalid BED format in file: " + bed_path +
                                 "\nOffending line: " + line);
    }

    if (start > end) {
        throw std::runtime_error("Start cannot be greater than end in BED file: " +
                                 bed_path + "\nOffending line: " + line);
    }

    Strand strand = strand_from_char(strand_str.empty() ? '.' : strand_str[0]);

    JunctionKey key{chrom, start, end, strand};
    junctions[key] = {name, score};
    return true;
}

std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash>
BedParser::parse(const std::string& bed_path) {
    std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash> junctions;

    if (is_gzipped(bed_path)) {
        // Use BGZF for gzipped files
        BGZF* fp = bgzf_open(bed_path.c_str(), "r");
        if (!fp) {
            throw std::runtime_error("Cannot open gzipped BED file: " + bed_path);
        }

        kstring_t str = {0, 0, nullptr};
        while (bgzf_getline(fp, '\n', &str) >= 0) {
            std::string line(str.s, str.l);
            process_bed_line(line, bed_path, junctions);
        }
        free(str.s);
        bgzf_close(fp);
    } else {
        // Use standard ifstream for plain text files
        std::ifstream file(bed_path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open BED file: " + bed_path);
        }

        std::string line;
        while (std::getline(file, line)) {
            process_bed_line(line, bed_path, junctions);
        }
    }

    return junctions;
}

JunctionMap BedParser::parse_multi(const std::vector<std::string>& bed_paths,
                                    int num_threads) {
    // Parse files in parallel
    std::vector<std::future<std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash>>> futures;

    for (const auto& path : bed_paths) {
        futures.push_back(std::async(
            num_threads > 1 ? std::launch::async : std::launch::deferred,
            [path]() { return parse(path); }
        ));
    }

    // Merge results
    JunctionMap merged;
    for (size_t i = 0; i < futures.size(); ++i) {
        auto result = futures[i].get();

        // Extract sample name from filename
        std::string sample_name;
        size_t last_slash = bed_paths[i].find_last_of("/\\");
        size_t last_dot = bed_paths[i].find_last_of('.');
        if (last_slash == std::string::npos) last_slash = 0;
        else last_slash++;
        if (last_dot == std::string::npos || last_dot < last_slash) {
            sample_name = bed_paths[i].substr(last_slash);
        } else {
            sample_name = bed_paths[i].substr(last_slash, last_dot - last_slash);
        }

        for (const auto& [key, value] : result) {
            auto& junction = merged[key];
            junction.add_sample(sample_name, value.second);
        }
    }

    return merged;
}

void BedParser::write(const std::string& output_path,
                      const JunctionMap& junctions,
                      const ScoringParams& scoring) {
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open output BED file: " + output_path);
    }

    // Sort junctions for consistent output
    std::vector<std::pair<JunctionKey, const JunctionData*>> sorted;
    for (const auto& [key, data] : junctions) {
        sorted.emplace_back(key, &data);
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    int junction_num = 1;
    for (const auto& [key, data] : sorted) {
        file << key.chrom << "\t"
             << key.start << "\t"
             << key.end << "\t"
             << "JUNC" << junction_num++ << "\t"
             << data->total_score << "\t"
             << to_char(key.strand) << "\n";
    }
}

} // namespace eastr
