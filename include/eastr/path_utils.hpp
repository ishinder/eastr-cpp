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

// Directory in which to place build artifacts (e.g. the auto-built bowtie2
// index) for the output `out_path`. A bare filename returns "" — the sentinel
// meaning "no explicit directory; use the default location". Crucially this
// never treats a single output *file* as a directory.
//   ""             -> ""
//   "filtered.bam" -> ""        (bare file: default / next to the reference)
//   "d/filt.bam"   -> "d"
//   "outdir"       -> "outdir"  (no extension: a directory)
std::string index_output_dir(const std::string& out_path);

// A real, never-empty directory in which to place a sidecar/temp file next to
// the output `out_path` (e.g. the temporary spurious-junction BED). Like
// index_output_dir but a bare filename resolves to "." (the current directory)
// rather than the empty sentinel, since a temp file needs a concrete location.
//   ""             -> "."
//   "filtered.bam" -> "."
//   "d/filt.bam"   -> "d"
//   "outdir"       -> "outdir"
std::string sidecar_dir(const std::string& out_path);

} // namespace eastr
