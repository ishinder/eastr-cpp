#include "eastr/bowtie2_runner.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <io.h>
#define popen _popen
#define pclose _pclose
#else
#include <unistd.h>
#endif

namespace eastr {

Bowtie2Runner::Bowtie2Runner(const std::string& bt2_index,
                             int num_threads,
                             int k_alignments)
    : bt2_index_(bt2_index),
      num_threads_(num_threads),
      k_alignments_(k_alignments) {
}

Bowtie2Runner::~Bowtie2Runner() = default;

std::string Bowtie2Runner::get_middle_seq(int target_len, const std::string& seq) {
    if (static_cast<int>(seq.size()) <= target_len) {
        return seq;
    }

    int middle_start = (static_cast<int>(seq.size()) - target_len) / 2;
    return seq.substr(middle_start, target_len);
}

std::string Bowtie2Runner::create_probe_fasta(
    const JunctionMap& introns,
    const SequenceCache& seqs,
    int overhang,
    int probe_length) {

    std::ostringstream fasta;
    int double_len = probe_length * 2;

    for (const auto& [key, data] : introns) {
        if (!data.hit) continue;

        auto rseq_it = seqs.find(data.jstart_region);
        auto qseq_it = seqs.find(data.jend_region);

        if (rseq_it == seqs.end() || qseq_it == seqs.end()) {
            continue;
        }

        const std::string& rseq = rseq_it->second;
        const std::string& qseq = qseq_it->second;
        const AlignmentHit& hit = *data.hit;

        // Create read name with junction coordinates
        std::string read_name = key.chrom + "," + std::to_string(key.start) + "," +
                                std::to_string(key.end) + "," + to_char(key.strand);

        // seq1: middle of aligned region from reference (5') flank
        if (hit.r_st < static_cast<int>(rseq.size()) &&
            hit.r_en <= static_cast<int>(rseq.size()) &&
            hit.r_st < hit.r_en) {
            std::string seq1 = get_middle_seq(double_len, rseq.substr(hit.r_st, hit.r_en - hit.r_st));
            fasta << ">" << read_name << ",seq1\n" << seq1 << "\n";
        }

        // seq2: middle of aligned region from query (3') flank
        if (hit.q_st < static_cast<int>(qseq.size()) &&
            hit.q_en <= static_cast<int>(qseq.size()) &&
            hit.q_st < hit.q_en) {
            std::string seq2 = get_middle_seq(double_len, qseq.substr(hit.q_st, hit.q_en - hit.q_st));
            fasta << ">" << read_name << ",seq2\n" << seq2 << "\n";
        }

        // seqh: hybrid sequence spanning the junction
        if (overhang >= probe_length &&
            static_cast<int>(rseq.size()) >= overhang &&
            static_cast<int>(qseq.size()) >= overhang + probe_length) {
            std::string seqh = rseq.substr(overhang - probe_length, probe_length) +
                               qseq.substr(overhang, probe_length);
            fasta << ">" << read_name << ",seqh\n" << seqh << "\n";
        }
    }

    return fasta.str();
}

void Bowtie2Runner::parse_sam_output(const std::string& sam_path,
                                     JunctionMap& introns) {
    std::ifstream sam_file(sam_path);
    if (!sam_file.is_open()) {
        throw std::runtime_error("Cannot open SAM file: " + sam_path);
    }

    std::string line;
    while (std::getline(sam_file, line)) {
        // Skip header lines
        if (line.empty() || line[0] == '@') continue;

        std::istringstream iss(line);
        std::string qname;
        int flag;

        if (!(iss >> qname >> flag)) continue;

        // Parse qname to get junction key and sequence type
        // Format: chrom,start,end,strand,seq_type
        std::vector<std::string> parts;
        std::istringstream qname_stream(qname);
        std::string part;
        while (std::getline(qname_stream, part, ',')) {
            parts.push_back(part);
        }

        if (parts.size() < 5) continue;

        JunctionKey key;
        key.chrom = parts[0];
        key.start = std::stoll(parts[1]);
        key.end = std::stoll(parts[2]);
        key.strand = strand_from_char(parts[3][0]);

        std::string seq_type = parts[4];

        // Find junction in map
        auto it = introns.find(key);
        if (it == introns.end()) continue;

        // Check if alignment is mapped (flag & 4 == unmapped)
        bool is_unmapped = (flag & 4) != 0;

        if (seq_type == "seqh") {
            // For hybrid sequence, count as 0 if unmapped
            if (is_unmapped) {
                it->second.seqh_count = 0;
            } else {
                it->second.seqh_count++;
            }
        } else if (seq_type == "seq1") {
            if (!is_unmapped) {
                it->second.seq1_count++;
            }
        } else if (seq_type == "seq2") {
            if (!is_unmapped) {
                it->second.seq2_count++;
            }
        }
    }
}

void Bowtie2Runner::align_probes(
    JunctionMap& introns_to_align,
    const SequenceCache& seqs,
    int overhang,
    int probe_length) {

    // Create temporary files
    char fasta_template[] = "/tmp/eastr_probes_XXXXXX";
    char sam_template[] = "/tmp/eastr_sam_XXXXXX";

    int fasta_fd = mkstemp(fasta_template);
    int sam_fd = mkstemp(sam_template);

    if (fasta_fd < 0 || sam_fd < 0) {
        throw std::runtime_error("Failed to create temporary files");
    }

    close(fasta_fd);
    close(sam_fd);

    std::string fasta_path = fasta_template;
    std::string sam_path = sam_template;

    // Write probe FASTA
    std::string fasta_content = create_probe_fasta(introns_to_align, seqs, overhang, probe_length);

    std::ofstream fasta_file(fasta_path);
    if (!fasta_file.is_open()) {
        throw std::runtime_error("Cannot write probe FASTA file");
    }
    fasta_file << fasta_content;
    fasta_file.close();

    // Build bowtie2 command
    // bowtie2 -p {p} --end-to-end -k {bt2_k} -D 20 -R 5 -L 20 -N 1 -i S,1,0.50
    std::ostringstream cmd;
    cmd << "bowtie2"
        << " -p " << num_threads_
        << " --end-to-end"
        << " -k " << k_alignments_
        << " -D 20 -R 5 -L 20 -N 1"
        << " -i S,1,0.50"
        << " -x " << bt2_index_
        << " -f " << fasta_path
        << " -S " << sam_path
        << " 2>/dev/null";

    // Run bowtie2
    int ret = system(cmd.str().c_str());
    if (ret != 0) {
        std::remove(fasta_path.c_str());
        std::remove(sam_path.c_str());
        throw std::runtime_error("bowtie2 command failed");
    }

    // Parse SAM output
    parse_sam_output(sam_path, introns_to_align);

    // Cleanup temporary files
    std::remove(fasta_path.c_str());
    std::remove(sam_path.c_str());
}

} // namespace eastr
