#pragma once

#include "junction.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace eastr {

class GtfParser {
public:
    // Internal transcript representation (public for helper functions)
    struct Transcript {
        std::string chrom;
        Strand strand;
        std::string gene_id;
        std::vector<std::pair<int64_t, int64_t>> exons;  // (start, end) pairs
    };

    // Extract splice junctions from GTF exon annotations
    // Supports both plain text and gzipped (.gz) GTF files
    static JunctionMap parse(const std::string& gtf_path);

    // Parse GTF attribute string (public for helper functions)
    static void parse_attributes(const std::string& attr_str,
                                 std::string& gene_id,
                                 std::string& transcript_id);
};

} // namespace eastr
