#pragma once

#include <string>
#include <vector>

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

enum class ListKind { Bam, Bed };

// Read a text file listing one input path per line (blank lines ignored).
// Throws std::runtime_error if the first path does not exist (the file is then
// probably neither a valid input nor a list), or if a later listed path is missing.
std::vector<std::string> read_path_list(const std::string& path, ListKind kind);

// Resolve a --bam/--bed argument into its input files: a single file, or the
// paths listed in it. The distinction is made by content (is_alignment_file /
// is_bed_file), not by name. `force_list` (--bam_list/--bed_list) skips
// detection and reads the argument as a list.
std::vector<std::string> resolve_input_files(const std::string& path, ListKind kind,
                                             bool force_list);

} // namespace eastr
