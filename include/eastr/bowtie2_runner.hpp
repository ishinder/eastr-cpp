#pragma once

#include "junction.hpp"

#include <string>
#include <unordered_map>

namespace eastr {

class Bowtie2Runner {
public:
    Bowtie2Runner(const std::string& bt2_index,
                  int num_threads,
                  int k_alignments);  // bt2_k + 1

    ~Bowtie2Runner();

    // Alignment counts for a junction
    struct AlignmentCounts {
        int seq1_count = 0;  // Alignments for ref overlap region
        int seq2_count = 0;  // Alignments for query overlap region
        int seqh_count = 0;  // Alignments for hybrid sequence
    };

    // Align probe sequences to reference
    // Updates junction data with alignment counts
    void align_probes(
        JunctionMap& introns_to_align,
        const SequenceCache& seqs,
        int overhang,
        int probe_length = 15);

private:
    std::string bt2_index_;
    int num_threads_;
    int k_alignments_;

    // Create probe FASTA content
    std::string create_probe_fasta(
        const JunctionMap& introns,
        const SequenceCache& seqs,
        int overhang,
        int probe_length);

    // Parse SAM output and update counts
    void parse_sam_output(
        const std::string& sam_path,
        JunctionMap& introns);

    // Get middle portion of sequence
    static std::string get_middle_seq(int target_len, const std::string& seq);
};

} // namespace eastr
