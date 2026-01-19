#include "eastr/gtf_parser.hpp"

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace eastr {

// Helper to check if file is gzipped
static bool is_gzipped(const std::string& path) {
    return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
}

void GtfParser::parse_attributes(const std::string& attr_str,
                                  std::string& gene_id,
                                  std::string& transcript_id) {
    gene_id.clear();
    transcript_id.clear();

    std::istringstream iss(attr_str);
    std::string token;

    while (std::getline(iss, token, ';')) {
        // Trim leading whitespace
        size_t start = token.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        token = token.substr(start);

        // Find attribute name and value
        size_t space_pos = token.find(' ');
        if (space_pos == std::string::npos) continue;

        std::string attr_name = token.substr(0, space_pos);
        std::string attr_value = token.substr(space_pos + 1);

        // Remove quotes from value
        if (!attr_value.empty() && attr_value.front() == '"') {
            attr_value = attr_value.substr(1);
        }
        if (!attr_value.empty() && attr_value.back() == '"') {
            attr_value.pop_back();
        }

        if (attr_name == "gene_id") {
            gene_id = attr_value;
        } else if (attr_name == "transcript_id") {
            transcript_id = attr_value;
        }
    }
}

// Process a single GTF line and add exon to transcripts map
static void process_gtf_line(const std::string& line, const std::string& gtf_path,
                             std::unordered_map<std::string, GtfParser::Transcript>& transcripts) {
    // Skip empty lines and comments
    if (line.empty() || line[0] == '#') return;

    std::istringstream iss(line);
    std::string chrom, source, feature;
    int64_t start, end;
    std::string score_str, strand_str, frame;

    // Parse GTF fields: chrom source feature start end score strand frame attributes
    if (!(iss >> chrom >> source >> feature >> start >> end >> score_str >> strand_str >> frame)) {
        return;  // Skip malformed lines
    }

    // Only process exon features
    if (feature != "exon") return;

    // Get the rest as attributes
    std::string attributes;
    std::getline(iss, attributes);

    // Parse attributes
    std::string gene_id, transcript_id;
    GtfParser::parse_attributes(attributes, gene_id, transcript_id);

    if (transcript_id.empty()) {
        throw std::runtime_error("Exon does not contain transcript_id in GTF file");
    }

    if (gene_id.empty()) {
        gene_id = "NA";
    }

    if (start > end) {
        throw std::runtime_error("Start cannot be greater than end in GTF file: " +
                                 gtf_path + "\nOffending line: " + line);
    }

    Strand strand = strand_from_char(strand_str.empty() ? '.' : strand_str[0]);

    // Add exon to transcript
    auto& trans = transcripts[transcript_id];
    if (trans.exons.empty()) {
        trans.chrom = chrom;
        trans.strand = strand;
        trans.gene_id = gene_id;
    }
    // GTF uses 1-based coordinates, convert to 0-based
    trans.exons.emplace_back(start - 1, end);
}

JunctionMap GtfParser::parse(const std::string& gtf_path) {
    std::unordered_map<std::string, Transcript> transcripts;

    if (is_gzipped(gtf_path)) {
        // Use BGZF for gzipped files
        BGZF* fp = bgzf_open(gtf_path.c_str(), "r");
        if (!fp) {
            throw std::runtime_error("Cannot open gzipped GTF file: " + gtf_path);
        }

        kstring_t str = {0, 0, nullptr};
        while (bgzf_getline(fp, '\n', &str) >= 0) {
            std::string line(str.s, str.l);
            process_gtf_line(line, gtf_path, transcripts);
        }
        free(str.s);
        bgzf_close(fp);
    } else {
        // Use standard ifstream for plain text files
        std::ifstream file(gtf_path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open GTF file: " + gtf_path);
        }

        std::string line;
        while (std::getline(file, line)) {
            process_gtf_line(line, gtf_path, transcripts);
        }
    }

    // Sort exons for each transcript and extract introns
    JunctionMap junctions;

    for (auto& [transcript_id, trans] : transcripts) {
        // Sort exons by start position
        std::sort(trans.exons.begin(), trans.exons.end());

        // Extract introns between consecutive exons
        for (size_t i = 1; i < trans.exons.size(); ++i) {
            // Intron: from end of previous exon to start of current exon
            int64_t intron_start = trans.exons[i - 1].second;  // Previous exon end (0-based)
            int64_t intron_end = trans.exons[i].first;          // Current exon start (0-based)

            JunctionKey key{trans.chrom, intron_start, intron_end, trans.strand};

            auto& junction = junctions[key];
            if (junction.gene_id.empty()) {
                junction.gene_id = trans.gene_id;
            }
            junction.transcript_ids.push_back(transcript_id);
        }
    }

    return junctions;
}

} // namespace eastr
