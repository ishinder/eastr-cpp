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

// True if htslib identifies `path` as BAM, CRAM or SAM (full format detection,
// stricter than is_alignment_file's magic-byte check: a gzip'd FASTA or a .csi
// index is not an alignment file). Used to pick files out of a directory.
bool is_alignment_format(const std::string& path);

enum class ListKind { Bam, Bed };

// The BAM/CRAM/SAM files directly inside `dir` (not recursive; hidden files and
// non-alignment files such as .bai indexes are skipped), sorted by path.
// Throws std::runtime_error if there are none.
std::vector<std::string> alignment_files_in_directory(const std::string& dir);

// Read a text file listing one input path per line (blank lines ignored).
// Throws std::runtime_error if the first path does not exist (the file is then
// probably neither a valid input nor a list), or if a later listed path is missing.
std::vector<std::string> read_path_list(const std::string& path, ListKind kind);

// Resolve a --bam/--bed argument into its input files: a single file, the
// paths listed in it, or (--bam only) the alignment files in a directory. The
// distinction is made by content (is_alignment_file / is_bed_file), not by name.
std::vector<std::string> resolve_input_files(const std::string& path, ListKind kind);

} // namespace eastr
