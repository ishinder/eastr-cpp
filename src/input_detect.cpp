#include "eastr/input_detect.hpp"

#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

namespace eastr {

namespace {

// Read up to `n` leading bytes of a file (fewer if the file is shorter).
std::string read_magic(const std::string& path, size_t n) {
    std::ifstream f(path, std::ios::binary);
    std::string buf(n, '\0');
    f.read(&buf[0], static_cast<std::streamsize>(n));
    buf.resize(static_cast<size_t>(f.gcount()));
    return buf;
}

bool has_gzip_magic(const std::string& magic) {
    return magic.size() >= 2 &&
           static_cast<unsigned char>(magic[0]) == 0x1f &&
           static_cast<unsigned char>(magic[1]) == 0x8b;
}

bool is_integer(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s) {
        if (!std::isdigit(static_cast<unsigned char>(c))) return false;
    }
    return true;
}

// fs::exists that never throws: a data line used as a "path" can exceed
// NAME_MAX/PATH_MAX, which makes the throwing overload raise ENAMETOOLONG.
bool path_exists(const std::string& p) {
    std::error_code ec;
    return fs::exists(p, ec);
}

// True if `line` looks like a SAM alignment record: >=11 tab-separated fields
// with integer FLAG (col 2), POS (col 4) and MAPQ (col 5).
bool looks_like_sam_record(const std::string& line) {
    std::vector<std::string> cols;
    std::istringstream iss(line);
    std::string col;
    while (std::getline(iss, col, '\t')) cols.push_back(col);
    return cols.size() >= 11 && is_integer(cols[1]) && is_integer(cols[3]) &&
           is_integer(cols[4]);
}

// Lower-cased extension of `path` (e.g. ".bam").
std::string lower_extension(const std::string& path) {
    fs::path p(path);
    std::string ext = p.extension().string();
    for (auto& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext;
}

// Error text for an input that is neither a recognized data file nor a list of
// existing paths. File content is checked first; the extension is used only as
// a hint about what the user intended (detection itself never relies on it).
std::string describe_unrecognized_input(const std::string& path,
                                        const std::string& first_line, ListKind kind) {
    const std::string q = "'" + path + "'";

    if (kind == ListKind::Bam) {
        const std::string ext = lower_extension(path);
        if (looks_like_sam_record(first_line) || ext == ".sam") {
            return q + " looks like a SAM file without a header (no @HD/@SQ lines). "
                   "EASTR needs the header to resolve reference sequence names.\n"
                   "Re-export it with its header, e.g. 'samtools view -h -o out.sam in.bam', "
                   "or pass the BAM/CRAM directly.";
        }
        if (ext == ".bam" || ext == ".cram") {
            const std::string fmt = (ext == ".bam") ? "BAM" : "CRAM";
            return q + " has a " + ext + " extension but is not a valid " + fmt +
                   " file (" + (ext == ".bam" ? "no BGZF" : "no CRAM") +
                   " signature at the start). It may be truncated, corrupt, or not " +
                   fmt + " at all.";
        }
    }

    std::string msg = q + " is not a " +
        (kind == ListKind::Bam ? "BAM/CRAM/SAM" : "BED") +
        " file, and its first line is not an existing file path.\n" +
        (kind == ListKind::Bam
             ? "--bam takes an alignment file, or a text file listing one alignment path per line."
             : "--bed takes a BED file, or a text file listing one BED path per line.");
    // Only echo the line when it plausibly is a (wrong) path, not a data record.
    if (first_line.find('\t') == std::string::npos && first_line.size() <= 512) {
        msg += "\nFirst line: " + first_line;
    }
    return msg;
}

} // namespace

bool is_alignment_file(const std::string& path) {
    std::string magic = read_magic(path, 4);

    // BAM is BGZF-compressed, which begins with the gzip magic bytes.
    if (has_gzip_magic(magic)) return true;
    // CRAM files begin with the literal "CRAM".
    if (magic.rfind("CRAM", 0) == 0) return true;
    // Uncompressed SAM begins with an '@' header line.
    if (!magic.empty() && magic[0] == '@') return true;

    // Otherwise assume a plain-text list of alignment paths.
    return false;
}

bool is_alignment_format(const std::string& path) {
    hFILE* hf = hopen(path.c_str(), "r");
    if (!hf) return false;
    htsFormat fmt;
    const bool ok = hts_detect_format2(hf, path.c_str(), &fmt) == 0;
    hclose_abruptly(hf);
    // Index files (.bai/.crai/.csi), FASTA, BED, VCF, plain gzip, ... are other
    // formats, so only real alignment files are selected.
    return ok && (fmt.format == bam || fmt.format == cram || fmt.format == sam);
}

bool is_bed_file(const std::string& path) {
    // A bgzipped BED file is a single data file, not a list.
    if (has_gzip_magic(read_magic(path, 2))) return true;

    std::ifstream f(path);
    if (!f.is_open()) {
        // Let downstream parsing surface a clear error; treat as single file.
        return true;
    }

    std::string line;
    while (std::getline(f, line)) {
        // Trim surrounding whitespace.
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty()) continue;
        if (line[0] == '#' || line.rfind("track", 0) == 0) continue;

        // First meaningful line found. A BED record has >=3 columns with
        // integer start/end (cols 2 and 3).
        std::istringstream iss(line);
        std::string chrom, start, end;
        if ((iss >> chrom >> start >> end) && is_integer(start) && is_integer(end)) {
            return true;
        }

        // Not a BED record: if it names an existing path, treat as a list file.
        if (path_exists(line)) {
            return false;
        }

        // Neither a BED record nor an existing path: treat as a (malformed) BED
        // file so the BED parser reports a precise format error.
        return true;
    }

    // Empty / comment-only file: treat as an empty BED file.
    return true;
}

std::vector<std::string> read_path_list(const std::string& path, ListKind kind) {
    std::ifstream list_file(path);
    if (!list_file.is_open()) {
        throw std::runtime_error("Cannot open file list: " + path);
    }

    std::vector<std::string> files;
    std::string line;
    while (std::getline(list_file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) continue;

        if (!path_exists(line)) {
            if (files.empty()) {
                // The first entry is not a path, so this is most likely not a
                // list at all. Diagnose what it probably is instead.
                throw std::runtime_error(describe_unrecognized_input(path, line, kind));
            }
            throw std::runtime_error("File not found: " + line + " (listed in " + path + ")");
        }
        files.push_back(line);
    }

    return files;
}

std::vector<std::string> alignment_files_in_directory(const std::string& dir) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        const std::string name = entry.path().filename().string();
        if (name.empty() || name[0] == '.') continue;  // hidden files (.DS_Store etc.)
        if (!entry.is_regular_file()) continue;        // not recursive; follows symlinks
        if (is_alignment_format(entry.path().string())) {
            files.push_back(entry.path().string());
        }
    }
    if (files.empty()) {
        throw std::runtime_error("'" + dir + "' is a directory with no BAM/CRAM/SAM files "
                                 "in it (subdirectories are not searched).");
    }
    // directory_iterator order is unspecified; sort so output order is stable.
    std::sort(files.begin(), files.end());
    return files;
}

std::vector<std::string> resolve_input_files(const std::string& path, ListKind kind,
                                             bool force_list) {
    if (!fs::exists(path)) {
        throw std::runtime_error("Input file not found: " + path);
    }

    if (fs::is_directory(path)) {
        if (kind != ListKind::Bam) {
            throw std::runtime_error("'" + path + "' is a directory. --bed takes a BED "
                                     "file, or a text file listing one BED path per line.");
        }
        if (force_list) {
            throw std::runtime_error("'" + path + "' is a directory, not a list file. "
                                     "Drop --bam_list to use the alignment files in it.");
        }
        return alignment_files_in_directory(path);
    }

    if (!force_list) {
        bool single = (kind == ListKind::Bam) ? is_alignment_file(path)
                                              : is_bed_file(path);
        if (single) {
            return {path};
        }
    } else if (kind == ListKind::Bam && is_alignment_file(path)) {
        // --bam_list on an actual BAM/CRAM/SAM: say so rather than misdiagnosing
        // its binary content as a broken list or corrupt file.
        throw std::runtime_error("'" + path + "' is an alignment file, not a list of "
                                 "paths. Drop --bam_list to process it directly.");
    }

    // A list of paths, either auto-detected or forced by --bam_list/--bed_list.
    return read_path_list(path, kind);
}

} // namespace eastr
