#pragma once

#include <string>

namespace eastr {

// Content-based input detection so input handling does not depend on a file's
// extension (Galaxy, for example, names every dataset "*.dat").

// True if `path` is a single alignment file (BAM/CRAM/SAM), detected by magic
// bytes: BGZF/gzip (1f 8b) -> BAM, "CRAM" -> CRAM, leading '@' -> SAM header.
// False means it should be treated as a text file listing alignment paths.
bool is_alignment_file(const std::string& path);

// True if `path` is a single BED file rather than a text file listing BED paths.
// A bgzipped file (1f 8b), or a first data line that parses as a BED record
// (>=3 whitespace columns with integer start/end), is a single BED file. If the
// first non-comment line is instead an existing filesystem path, it is a list.
bool is_bed_file(const std::string& path);

} // namespace eastr
