#include "eastr/alignment_io.hpp"

#include <htslib/hts.h>

namespace eastr {

namespace {
// Written once in main() before any worker threads start; read-only afterwards.
std::string g_reference_fasta;
} // namespace

void set_alignment_reference(const std::string& fasta_path) {
    g_reference_fasta = fasta_path;
}

samFile* open_alignment(const std::string& path, int threads) {
    samFile* fp = sam_open(path.c_str(), "r");
    if (!fp) return nullptr;

    if (fp->format.format == cram) {
        if (!g_reference_fasta.empty() &&
            hts_set_fai_filename(fp, g_reference_fasta.c_str()) != 0) {
            sam_close(fp);
            return nullptr;
        }
        // Decode on htslib pool threads (large stacks); see alignment_io.hpp.
        if (hts_set_threads(fp, threads > 1 ? threads : 1) != 0) {
            sam_close(fp);
            return nullptr;
        }
    } else if (threads > 1) {
        hts_set_threads(fp, threads);
    }
    return fp;
}

} // namespace eastr
