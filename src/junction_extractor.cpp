#include "eastr/junction_extractor.hpp"

#include <htslib/sam.h>

#include <algorithm>
#include <fstream>
#include <future>
#include <iomanip>
#include <stdexcept>

namespace eastr {

struct JunctionExtractor::Impl {
    std::string bam_path;
    std::string sample_name;
    samFile* bam_file = nullptr;
    sam_hdr_t* header = nullptr;
    std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash> junctions;

    ~Impl() {
        if (header) {
            sam_hdr_destroy(header);
        }
        if (bam_file) {
            sam_close(bam_file);
        }
    }
};

JunctionExtractor::JunctionExtractor(const std::string& bam_path)
    : pImpl_(std::make_unique<Impl>()) {
    pImpl_->bam_path = bam_path;

    // Extract sample name from filename
    size_t last_slash = bam_path.find_last_of("/\\");
    size_t last_dot = bam_path.find_last_of('.');
    if (last_slash == std::string::npos) last_slash = 0;
    else last_slash++;
    if (last_dot == std::string::npos || last_dot < last_slash) {
        pImpl_->sample_name = bam_path.substr(last_slash);
    } else {
        pImpl_->sample_name = bam_path.substr(last_slash, last_dot - last_slash);
    }

    // Open BAM file
    pImpl_->bam_file = sam_open(bam_path.c_str(), "r");
    if (!pImpl_->bam_file) {
        throw std::runtime_error("Cannot open BAM file: " + bam_path);
    }

    // Read header
    pImpl_->header = sam_hdr_read(pImpl_->bam_file);
    if (!pImpl_->header) {
        throw std::runtime_error("Cannot read BAM header: " + bam_path);
    }
}

JunctionExtractor::~JunctionExtractor() = default;

JunctionExtractor::JunctionExtractor(JunctionExtractor&&) noexcept = default;
JunctionExtractor& JunctionExtractor::operator=(JunctionExtractor&&) noexcept = default;

const std::string& JunctionExtractor::sample_name() const {
    return pImpl_->sample_name;
}

std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash>
JunctionExtractor::extract() {
    std::unordered_map<JunctionKey, int, JunctionKeyHash> junction_counts;

    bam1_t* aln = bam_init1();
    if (!aln) {
        throw std::runtime_error("Failed to allocate BAM alignment");
    }

    while (sam_read1(pImpl_->bam_file, pImpl_->header, aln) >= 0) {
        // Skip unmapped reads
        if (aln->core.flag & BAM_FUNMAP) continue;

        // Get chromosome name
        const char* chrom = sam_hdr_tid2name(pImpl_->header, aln->core.tid);
        if (!chrom) continue;

        // Get splice strand from XS tag (set by aligners like HISAT2)
        // This is the transcription strand, not the alignment strand
        Strand strand = Strand::Unknown;
        uint8_t* xs_tag = bam_aux_get(aln, "XS");
        if (xs_tag) {
            char xs_char = bam_aux2A(xs_tag);
            if (xs_char == '+') strand = Strand::Plus;
            else if (xs_char == '-') strand = Strand::Minus;
        }

        // Get YC tag for duplicate count (used by some pipelines like StringTie)
        // Default to 1 if not present
        int yc_count = 1;
        uint8_t* yc_tag = bam_aux_get(aln, "YC");
        if (yc_tag) {
            yc_count = bam_aux2i(yc_tag);
            if (yc_count < 1) yc_count = 1;  // Ensure at least 1
        }

        // Parse CIGAR to find splice junctions (N operations)
        uint32_t* cigar = bam_get_cigar(aln);
        int64_t ref_pos = aln->core.pos;  // 0-based position

        for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);

            if (op == BAM_CREF_SKIP) {  // N operation - splice junction
                // Junction spans from ref_pos to ref_pos + len
                JunctionKey key{chrom, ref_pos, ref_pos + len, strand};
                junction_counts[key] += yc_count;  // Use YC tag count
            }

            // Update reference position based on operation
            if (bam_cigar_type(op) & 2) {  // Consumes reference
                ref_pos += len;
            }
        }
    }

    bam_destroy1(aln);

    // Convert counts to result format
    std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash> result;
    for (const auto& [key, count] : junction_counts) {
        result[key] = {pImpl_->sample_name, count};
    }

    pImpl_->junctions = result;
    return result;
}

void JunctionExtractor::write_bed(const std::string& output_path) {
    if (pImpl_->junctions.empty()) {
        extract();
    }

    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open output BED file: " + output_path);
    }

    // No track header line - matches reference behavior

    // Sort junctions for consistent output
    std::vector<std::pair<JunctionKey, std::pair<std::string, int>>> sorted(
        pImpl_->junctions.begin(), pImpl_->junctions.end());
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    int junction_num = 1;
    for (const auto& [key, value] : sorted) {
        // Use 8-digit zero-padded junction names to match reference format
        file << key.chrom << "\t"
             << key.start << "\t"
             << key.end << "\t"
             << "JUNC" << std::setfill('0') << std::setw(8) << junction_num++ << "\t"
             << value.second << "\t"
             << to_char(key.strand) << "\n";
    }
}

JunctionMap extract_junctions_multi_bam(
    const std::vector<std::string>& bam_paths,
    const std::vector<std::string>& output_paths,
    int num_threads) {

    // Define extraction task
    auto extract_task = [](const std::string& bam_path, const std::string& output_path) {
        JunctionExtractor extractor(bam_path);
        auto junctions = extractor.extract();
        if (!output_path.empty()) {
            extractor.write_bed(output_path);
        }
        return std::make_pair(extractor.sample_name(), junctions);
    };

    // Launch extraction tasks
    std::vector<std::future<std::pair<std::string,
        std::unordered_map<JunctionKey, std::pair<std::string, int>, JunctionKeyHash>>>> futures;

    for (size_t i = 0; i < bam_paths.size(); ++i) {
        std::string output = (i < output_paths.size()) ? output_paths[i] : "";
        futures.push_back(std::async(
            num_threads > 1 ? std::launch::async : std::launch::deferred,
            extract_task, bam_paths[i], output
        ));
    }

    // Merge results
    JunctionMap merged;
    for (auto& future : futures) {
        auto [sample_name, junctions] = future.get();

        for (const auto& [key, value] : junctions) {
            auto& junction = merged[key];
            junction.add_sample(value.first, value.second);
        }
    }

    return merged;
}

} // namespace eastr
