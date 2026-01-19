#include "eastr/junction_accumulator.hpp"

#include <algorithm>

namespace eastr {

JunctionAccumulator::JunctionAccumulator(size_t num_files,
                                         const std::vector<std::string>& sample_names)
    : sample_names_(sample_names) {
    file_states_.resize(num_files);
}

void JunctionAccumulator::process_alignment(const MultiBamReader::AlignmentRecord& rec,
                                            sam_hdr_t* header) {
    bam1_t* aln = rec.aln;

    // Skip unmapped reads
    if (aln->core.flag & BAM_FUNMAP) {
        return;
    }

    // Update file position tracking
    const char* chrom = sam_hdr_tid2name(header, aln->core.tid);
    if (!chrom) return;

    std::string chrom_str(chrom);
    update_file_position(rec.file_index, chrom_str, aln->core.pos);

    // Check if this is a spliced alignment (has N in CIGAR)
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

    // Extract junctions from this alignment
    auto junctions = extract_junctions_from_alignment(rec, header);

    // Add to active junctions
    for (auto& [key, score] : junctions) {
        auto& junction = active_junctions_[key];

        // Initialize if new
        if (junction.total_score == 0) {
            total_junctions_++;
        }

        // Add this sample's contribution
        junction.add_sample(rec.sample_name, score);
    }
}

std::vector<std::pair<JunctionKey, int>>
JunctionAccumulator::extract_junctions_from_alignment(
    const MultiBamReader::AlignmentRecord& rec,
    sam_hdr_t* header) {

    std::vector<std::pair<JunctionKey, int>> result;

    bam1_t* aln = rec.aln;
    const char* chrom = sam_hdr_tid2name(header, aln->core.tid);
    if (!chrom) return result;

    std::string chrom_str(chrom);

    // Get splice strand from XS tag (set by aligners like HISAT2)
    // This is the transcription strand, not the alignment strand
    // Match original junction_extractor behavior - only use XS tag
    Strand strand = Strand::Unknown;
    uint8_t* xs_tag = bam_aux_get(aln, "XS");
    if (xs_tag) {
        char xs_val = bam_aux2A(xs_tag);
        if (xs_val == '+') strand = Strand::Plus;
        else if (xs_val == '-') strand = Strand::Minus;
    }

    // Parse CIGAR to find splice junctions
    uint32_t* cigar = bam_get_cigar(aln);
    int64_t ref_pos = aln->core.pos;  // 0-based

    for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CREF_SKIP) {
            // This is a splice junction (intron)
            // Junction: [ref_pos, ref_pos + len) in 0-based half-open
            JunctionKey key;
            key.chrom = chrom_str;
            key.start = ref_pos;
            key.end = ref_pos + len;
            key.strand = strand;

            result.emplace_back(key, 1);  // Score = 1 for each supporting read
        }

        // Update reference position for operations that consume reference
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
            op == BAM_CEQUAL || op == BAM_CDIFF) {
            ref_pos += len;
        }
    }

    return result;
}

std::vector<std::pair<JunctionKey, JunctionData>>
JunctionAccumulator::get_completed_junctions() {
    std::vector<std::pair<JunctionKey, JunctionData>> completed;

    // Find junctions where all files have advanced past their end position
    auto it = active_junctions_.begin();
    while (it != active_junctions_.end()) {
        const auto& key = it->first;

        if (is_position_complete(key.chrom, key.end)) {
            completed.emplace_back(key, std::move(it->second));
            it = active_junctions_.erase(it);
        } else {
            ++it;
        }
    }

    return completed;
}

std::vector<std::pair<JunctionKey, JunctionData>>
JunctionAccumulator::flush() {
    std::vector<std::pair<JunctionKey, JunctionData>> result;
    result.reserve(active_junctions_.size());

    for (auto& [key, data] : active_junctions_) {
        result.emplace_back(key, std::move(data));
    }

    active_junctions_.clear();
    return result;
}

void JunctionAccumulator::notify_chromosome_change(int file_idx, int new_tid) {
    if (file_idx >= 0 && static_cast<size_t>(file_idx) < file_states_.size()) {
        if (new_tid < 0) {
            file_states_[file_idx].exhausted = true;
        }
        file_states_[file_idx].min_pos = -1;  // Reset position for new chromosome
    }
}

bool JunctionAccumulator::is_position_complete(const std::string& chrom, int64_t end_pos) const {
    for (const auto& state : file_states_) {
        if (state.exhausted) {
            continue;  // Exhausted files don't block completion
        }

        // If file is on a different (later) chromosome, this position is complete
        if (state.current_chrom != chrom) {
            // Check if file has moved to a later chromosome
            // For simplicity, assume files process chromosomes in same order
            // If file is still on earlier chrom or hasn't seen this chrom, not complete
            if (state.current_chrom.empty() || state.current_chrom < chrom) {
                return false;
            }
            // File is on later chromosome, so this position is complete for this file
            continue;
        }

        // Same chromosome - check if position has advanced past end_pos
        if (state.min_pos < end_pos) {
            return false;
        }
    }

    return true;
}

void JunctionAccumulator::update_file_position(int file_idx, const std::string& chrom, int64_t pos) {
    if (file_idx < 0 || static_cast<size_t>(file_idx) >= file_states_.size()) {
        return;
    }

    auto& state = file_states_[file_idx];

    if (state.current_chrom != chrom) {
        state.current_chrom = chrom;
        state.min_pos = pos;
    } else {
        // Only update if position advanced (should always be true for sorted BAMs)
        if (pos > state.min_pos) {
            state.min_pos = pos;
        }
    }

    // Update global current chromosome tracking
    if (current_chrom_.empty() || current_chrom_ != chrom) {
        // Check if all files are now on this chromosome or later
        bool all_advanced = true;
        for (const auto& s : file_states_) {
            if (!s.exhausted && (s.current_chrom.empty() || s.current_chrom < chrom)) {
                all_advanced = false;
                break;
            }
        }
        if (all_advanced) {
            current_chrom_ = chrom;
        }
    }
}

} // namespace eastr
