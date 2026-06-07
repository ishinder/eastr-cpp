#include "eastr/input_detect.hpp"

#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

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
        if (fs::exists(line)) {
            return false;
        }

        // Neither a BED record nor an existing path: treat as a (malformed) BED
        // file so the BED parser reports a precise format error.
        return true;
    }

    // Empty / comment-only file: treat as an empty BED file.
    return true;
}

} // namespace eastr
