#pragma once

#include "junction.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace eastr {

class BedParser {
public:
    // Parse single BED file
    // Returns map of junction key -> (name, score)
    static std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash>
    parse(const std::string& bed_path);

    // Parse multiple BED files and merge into JunctionMap
    static JunctionMap parse_multi(
        const std::vector<std::string>& bed_paths,
        int num_threads = 1);

    // Write junctions to BED file
    static void write(const std::string& output_path,
                      const JunctionMap& junctions,
                      const ScoringParams& scoring);
};

} // namespace eastr
