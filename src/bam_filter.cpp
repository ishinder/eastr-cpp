#include "eastr/bam_filter.hpp"
#include "eastr/bed_parser.hpp"

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <future>
#include <iostream>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstring>

namespace eastr {

// Key for tracking alignments: (read_name, ref_name, pos, mate_pos)
using AlnKey = std::tuple<std::string, std::string, int64_t, int64_t>;

// Stored alignment record for exact matching
struct StoredAlignment {
    std::vector<uint32_t> cigar;
    int pair_order;  // 1 or 2
};

struct BamFilter::Impl {
    // NH tag tracking: key = "readname;pairorder", value = count to subtract
    std::unordered_map<std::string, int> nh_adjustments;

    // Removed alignments for exact matching
    std::map<AlnKey, std::vector<StoredAlignment>> removed_alignments;
};

BamFilter::BamFilter(const Options& options)
    : options_(options), pImpl_(std::make_unique<Impl>()) {
}

BamFilter::~BamFilter() = default;

std::unordered_set<JunctionKey, JunctionKeyHash>
BamFilter::load_spurious_junctions(const std::string& bed_path) {
    auto junctions = BedParser::parse(bed_path);

    std::unordered_set<JunctionKey, JunctionKeyHash> spurious;
    for (const auto& [key, _] : junctions) {
        spurious.insert(key);
    }

    return spurious;
}

// Check if two alignments have identical CIGAR
static bool check_identical_cigar(const bam1_t* rec1, const std::vector<uint32_t>& cigar2) {
    if (rec1->core.n_cigar != cigar2.size()) {
        return false;
    }
    return memcmp(bam_get_cigar(rec1), cigar2.data(),
                  rec1->core.n_cigar * sizeof(uint32_t)) == 0;
}

// Get pair order (1 or 2) from BAM record
static int get_pair_order(const bam1_t* aln) {
    if (aln->core.flag & BAM_FREAD1) return 1;
    if (aln->core.flag & BAM_FREAD2) return 2;
    return 0;
}

// Create NH adjustment key
static std::string make_nh_key(const char* qname, int pair_order) {
    return std::string(qname) + ";" + std::to_string(pair_order);
}

BamFilter::FilterStats BamFilter::filter(
    const std::string& input_bam,
    const std::string& spurious_bed,
    const std::string& output_bam) {

    FilterStats stats;
    pImpl_->nh_adjustments.clear();
    pImpl_->removed_alignments.clear();

    // Load spurious junctions
    auto spurious = load_spurious_junctions(spurious_bed);

    // Open input BAM
    samFile* in = sam_open(input_bam.c_str(), "r");
    if (!in) {
        throw std::runtime_error("Cannot open input BAM: " + input_bam);
    }

    // Set I/O threads for BGZF decompression
    if (options_.io_threads > 1) {
        hts_set_threads(in, options_.io_threads);
    }

    // Read header
    sam_hdr_t* header = sam_hdr_read(in);
    if (!header) {
        sam_close(in);
        throw std::runtime_error("Cannot read BAM header: " + input_bam);
    }

    bam1_t* aln = bam_init1();
    int64_t num_unmapped = 0;
    int64_t num_removed_spliced = 0;

    // ============ FIRST PASS: Identify alignments to remove ============
    while (sam_read1(in, header, aln) >= 0) {
        stats.total_alignments++;

        // Skip unmapped reads
        if (aln->core.flag & BAM_FUNMAP) {
            num_unmapped++;
            continue;
        }

        // Get chromosome name
        const char* chrom = sam_hdr_tid2name(header, aln->core.tid);
        if (!chrom) continue;

        // Check for splice junctions (N operations in CIGAR)
        uint32_t* cigar = bam_get_cigar(aln);
        int64_t ref_pos = aln->core.pos;
        bool has_spurious = false;
        bool is_spliced = false;

        // Get strand from XS tag or FLAG
        char strand_char = '.';
        uint8_t* xs_tag = bam_aux_get(aln, "XS");
        if (xs_tag) {
            strand_char = bam_aux2A(xs_tag);
        }
        Strand strand = (strand_char == '-') ? Strand::Minus :
                       (strand_char == '+') ? Strand::Plus : Strand::Unknown;

        for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);

            if (op == BAM_CREF_SKIP) {  // N operation - splice junction
                is_spliced = true;
                stats.total_spliced++;

                // Junction coordinates: [ref_pos, ref_pos + len)
                // Note: vacuum uses [exon_end, next_exon_start-1]
                JunctionKey key{chrom, ref_pos, ref_pos + len, strand};
                if (spurious.find(key) != spurious.end()) {
                    has_spurious = true;
                }
            }

            if (bam_cigar_type(op) & 2) {  // Consumes reference
                ref_pos += len;
            }
        }

        if (has_spurious) {
            num_removed_spliced++;
            const char* qname = bam_get_qname(aln);
            int pair_order = get_pair_order(aln);

            // Track NH tag adjustment
            std::string nh_key = make_nh_key(qname, pair_order);
            pImpl_->nh_adjustments[nh_key]++;

            // Store alignment for exact matching
            AlnKey aln_key = std::make_tuple(
                std::string(qname),
                std::string(chrom),
                static_cast<int64_t>(aln->core.pos),
                static_cast<int64_t>(aln->core.mpos)
            );

            StoredAlignment stored;
            stored.cigar.assign(cigar, cigar + aln->core.n_cigar);
            stored.pair_order = pair_order;
            pImpl_->removed_alignments[aln_key].push_back(stored);
        }
    }

    if (options_.verbose) {
        std::cerr << "First pass complete:\n"
                  << "  Total alignments: " << stats.total_alignments << "\n"
                  << "  Unmapped: " << num_unmapped << "\n"
                  << "  Spliced: " << stats.total_spliced << "\n"
                  << "  Flagged for removal: " << num_removed_spliced << "\n";
    }

    // ============ SECOND PASS: Write filtered output ============
    // Rewind the BAM file - always reopen to ensure clean state
    sam_close(in);
    in = sam_open(input_bam.c_str(), "r");
    if (!in) {
        bam_destroy1(aln);
        sam_hdr_destroy(header);
        throw std::runtime_error("Cannot reopen input BAM: " + input_bam);
    }

    // Set I/O threads for BGZF decompression
    if (options_.io_threads > 1) {
        hts_set_threads(in, options_.io_threads);
    }

    sam_hdr_destroy(header);
    header = sam_hdr_read(in);
    if (!header) {
        sam_close(in);
        bam_destroy1(aln);
        throw std::runtime_error("Cannot re-read BAM header: " + input_bam);
    }

    // Open output BAM
    samFile* out = sam_open(output_bam.c_str(), "wb");
    if (!out) {
        bam_destroy1(aln);
        sam_hdr_destroy(header);
        sam_close(in);
        throw std::runtime_error("Cannot open output BAM: " + output_bam);
    }

    // Set I/O threads for BGZF compression
    if (options_.io_threads > 1) {
        hts_set_threads(out, options_.io_threads);
    }

    // Write header
    if (sam_hdr_write(out, header) < 0) {
        sam_close(out);
        bam_destroy1(aln);
        sam_hdr_destroy(header);
        sam_close(in);
        throw std::runtime_error("Cannot write BAM header: " + output_bam);
    }

    // Optional: open removed alignments BAM
    samFile* removed_out = nullptr;
    if (options_.write_removed && !options_.removed_output_path.empty()) {
        removed_out = sam_open(options_.removed_output_path.c_str(), "wb");
        if (removed_out) {
            // Set I/O threads for BGZF compression
            if (options_.io_threads > 1) {
                hts_set_threads(removed_out, options_.io_threads);
            }
            sam_hdr_write(removed_out, header);
        }
    }

    // Track mates that have been unpaired
    std::map<AlnKey, int> mates_unpaired;
    int64_t num_mates_affected = 0;
    int64_t num_output = 0;

    while (sam_read1(in, header, aln) >= 0) {
        // Skip unmapped reads - don't write to output (matches Python behavior)
        if (aln->core.flag & BAM_FUNMAP) {
            continue;
        }

        const char* qname = bam_get_qname(aln);
        const char* chrom = (aln->core.tid >= 0) ?
                            sam_hdr_tid2name(header, aln->core.tid) : "*";

        // Check if this alignment should be removed
        AlnKey aln_key = std::make_tuple(
            std::string(qname),
            std::string(chrom),
            static_cast<int64_t>(aln->core.pos),
            static_cast<int64_t>(aln->core.mpos)
        );

        auto it = pImpl_->removed_alignments.find(aln_key);
        if (it != pImpl_->removed_alignments.end()) {
            // Check for exact CIGAR match
            bool found = false;
            for (const auto& stored : it->second) {
                if (check_identical_cigar(aln, stored.cigar)) {
                    found = true;
                    stats.removed_alignments++;
                    if (removed_out) {
                        sam_write1(removed_out, header, aln);
                    }
                    break;
                }
            }
            if (found) {
                continue;  // Skip this alignment
            }
        }

        // Check if this is a mate of a removed alignment
        if (aln->core.flag & BAM_FPAIRED) {
            AlnKey mate_key = std::make_tuple(
                std::string(qname),
                std::string(chrom),
                static_cast<int64_t>(aln->core.mpos),
                static_cast<int64_t>(aln->core.pos)
            );

            auto mate_it = pImpl_->removed_alignments.find(mate_key);
            if (mate_it != pImpl_->removed_alignments.end()) {
                int num_removed = mate_it->second.size();
                int& num_seen = mates_unpaired[mate_key];

                bool should_update = (num_seen < num_removed);
                if (should_update) {
                    num_seen++;
                    num_mates_affected++;
                }

                if (options_.remove_mate && should_update) {
                    // Remove the mate entirely
                    stats.removed_alignments++;
                    if (removed_out) {
                        sam_write1(removed_out, header, aln);
                    }
                    continue;
                }

                // Update NH tag if needed
                int pair_order = get_pair_order(aln);
                std::string nh_key = make_nh_key(qname, pair_order);
                auto nh_it = pImpl_->nh_adjustments.find(nh_key);
                if (nh_it != pImpl_->nh_adjustments.end()) {
                    uint8_t* nh_tag = bam_aux_get(aln, "NH");
                    if (nh_tag) {
                        int old_nh = bam_aux2i(nh_tag);
                        int new_nh = old_nh - nh_it->second;
                        if (new_nh > 0) {
                            bam_aux_update_int(aln, "NH", new_nh);
                        }
                    }
                }

                // Unpair the mate (update FLAG, TLEN, MPOS)
                if (should_update) {
                    aln->core.flag &= ~(BAM_FPAIRED | BAM_FPROPER_PAIR | BAM_FMUNMAP);
                    aln->core.isize = 0;  // Set template length to zero
                    aln->core.mpos = aln->core.pos;  // Set mate pos to own pos
                }
            }
        }

        // Write to output
        sam_write1(out, header, aln);
        num_output++;
    }

    // Cleanup
    bam_destroy1(aln);
    sam_hdr_destroy(header);
    sam_close(in);
    sam_close(out);

    if (removed_out) {
        sam_close(removed_out);
    }

    if (options_.verbose) {
        std::cerr << "Second pass complete:\n"
                  << "  Removed alignments: " << stats.removed_alignments << "\n"
                  << "  Mates affected: " << num_mates_affected << "\n"
                  << "  Alignments written: " << num_output << "\n";
    }

    return stats;
}

void filter_multi_bam(
    const std::vector<std::string>& input_bams,
    const std::string& spurious_bed,
    const std::vector<std::string>& output_bams,
    int num_threads,
    const BamFilter::Options& options) {

    // Calculate I/O threads per file to utilize available cores
    int io_threads_per_file = std::max(1, num_threads / static_cast<int>(input_bams.size()));

    // Define filter task - all samples use the same spurious BED file
    auto filter_task = [&options, &spurious_bed, io_threads_per_file](
                                   const std::string& input,
                                   const std::string& output) {
        BamFilter::Options task_options = options;
        task_options.io_threads = io_threads_per_file;

        // Set removed output path if write_removed is enabled
        if (task_options.write_removed) {
            // Generate path: output.bam -> output_removed_alignments.bam
            std::string removed_path = output;
            size_t dot_pos = removed_path.rfind(".bam");
            if (dot_pos != std::string::npos) {
                removed_path = removed_path.substr(0, dot_pos) + "_removed_alignments.bam";
            } else {
                removed_path += "_removed_alignments.bam";
            }
            task_options.removed_output_path = removed_path;
        }

        BamFilter filter(task_options);
        return filter.filter(input, spurious_bed, output);
    };

    // Launch filter tasks
    std::vector<std::future<BamFilter::FilterStats>> futures;

    for (size_t i = 0; i < input_bams.size(); ++i) {
        futures.push_back(std::async(
            num_threads > 1 ? std::launch::async : std::launch::deferred,
            filter_task, input_bams[i], output_bams[i]
        ));
    }

    // Wait for all tasks
    for (auto& future : futures) {
        future.get();
    }
}

} // namespace eastr
