#pragma once

#include <string>

namespace eastr {

// Decide whether `path` denotes a file (vs a directory) for output purposes.
// An existing directory, or a path with a trailing separator, is a directory.
// Otherwise it is treated as a file when it has an extension after the last
// path separator (e.g. "out.bed", "dir/out.bed"); a bare name without an
// extension (e.g. "results", "out_dir") is treated as a directory.
bool is_file_path(const std::string& path);

// Filename component of `path` with its directory and extension removed
// (e.g. "/a/b/sample.bam" -> "sample").
std::string path_basename(const std::string& path);

// Directory component of `path`, or "." if `path` has no separator.
std::string path_parent_dir(const std::string& path);

} // namespace eastr
