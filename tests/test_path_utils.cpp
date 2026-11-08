#include <catch2/catch_test_macros.hpp>

#include "eastr/path_utils.hpp"
#include "eastr/output_writer.hpp"

using namespace eastr;

TEST_CASE("is_file_path distinguishes files from directories", "[path_utils]") {
    // Bare filename with extension and no separator -> file (regression: this used
    // to be misclassified as a directory because of the npos comparison bug).
    REQUIRE(is_file_path("out.bed"));
    REQUIRE(is_file_path("out.bam"));

    // File with a directory prefix -> file.
    REQUIRE(is_file_path("dir/out.bed"));
    REQUIRE(is_file_path("a/b/c/out.original.bed"));

    // No extension -> directory.
    REQUIRE_FALSE(is_file_path("results"));
    REQUIRE_FALSE(is_file_path("out/results"));

    // Trailing separator -> directory.
    REQUIRE_FALSE(is_file_path("results/"));

    // Empty -> not a file.
    REQUIRE_FALSE(is_file_path(""));
}

TEST_CASE("path_basename strips directory and extension", "[path_utils]") {
    REQUIRE(path_basename("/a/b/sample.bam") == "sample");
    REQUIRE(path_basename("sample.bam") == "sample");
    REQUIRE(path_basename("sample") == "sample");
    REQUIRE(path_basename("dir/sample") == "sample");
}

TEST_CASE("path_parent_dir returns the directory component", "[path_utils]") {
    REQUIRE(path_parent_dir("/a/b/sample.bam") == "/a/b");
    REQUIRE(path_parent_dir("sample.bam") == ".");
}

// NOTE: index_output_dir/sidecar_dir call is_file_path, which inspects the
// filesystem (fs::is_directory). These cases use paths that do NOT exist on
// disk; do not create fixtures named "filtered.bam"/"outdir" or results flip.

TEST_CASE("index_output_dir: bare file -> default (empty), never a directory",
          "[path_utils]") {
    // Regression (2.1.1): a bare "filtered.bam" must NOT be used as a directory.
    // Empty result => build the index in the default location (next to reference).
    REQUIRE(index_output_dir("") == "");
    REQUIRE(index_output_dir("filtered.bam") == "");
    REQUIRE(index_output_dir("out/filtered.bam") == "out");
    REQUIRE(index_output_dir("a/b/c/f.bam") == "a/b/c");
    REQUIRE(index_output_dir("outdir") == "outdir");  // no extension => a directory
}

TEST_CASE("sidecar_dir: bare file -> '.', never a directory", "[path_utils]") {
    // Regression (2.1.1): a bare "filtered.bam" must resolve to "." (cwd), not
    // create a directory named "filtered.bam/".
    REQUIRE(sidecar_dir("") == ".");
    REQUIRE(sidecar_dir("filtered.bam") == ".");
    REQUIRE(sidecar_dir("out/filtered.bam") == "out");
    REQUIRE(sidecar_dir("a/b/filt.bam") == "a/b");
    REQUIRE(sidecar_dir("outdir") == "outdir");
}

TEST_CASE("generate_junction_output_paths: single input + file path", "[output_writer]") {
    // Single input, bare filename -> use the path as-is (issue #3 core).
    REQUIRE(OutputWriter::generate_junction_output_paths(
        {"sample.bam"}, "out.bed", "_original_junctions") ==
        std::vector<std::string>{"out.bed"});
}

TEST_CASE("generate_junction_output_paths: single input + directory", "[output_writer]") {
    REQUIRE(OutputWriter::generate_junction_output_paths(
        {"/data/sample.bam"}, "outdir", "_original_junctions") ==
        std::vector<std::string>{"outdir/sample_original_junctions.bed"});
}

TEST_CASE("generate_junction_output_paths: multi input falls back to parent dir",
          "[output_writer]") {
    REQUIRE(OutputWriter::generate_junction_output_paths(
        {"a.bam", "b.bam"}, "results/out.bed", "_original_junctions") ==
        std::vector<std::string>{"results/a_original_junctions.bed",
                                 "results/b_original_junctions.bed"});
}

TEST_CASE("generate_bam_output_paths: single input + file path", "[output_writer]") {
    // Single BAM, bare filename -> write that file directly (issue #4).
    REQUIRE(OutputWriter::generate_bam_output_paths(
        {"sample.bam"}, "out.bam", "_EASTR_filtered") ==
        std::vector<std::string>{"out.bam"});
}

TEST_CASE("generate_bam_output_paths: single input + directory", "[output_writer]") {
    // Single BAM + a directory path -> one per-sample file inside the directory.
    REQUIRE(OutputWriter::generate_bam_output_paths(
        {"/d/sample.bam"}, "outdir", "_EASTR_filtered") ==
        std::vector<std::string>{"outdir/sample_EASTR_filtered.bam"});
}

TEST_CASE("generate_bam_output_paths: multi input uses directory", "[output_writer]") {
    REQUIRE(OutputWriter::generate_bam_output_paths(
        {"a.bam", "b.bam"}, "outdir", "_EASTR_filtered") ==
        std::vector<std::string>{"outdir/a_EASTR_filtered.bam",
                                 "outdir/b_EASTR_filtered.bam"});
}

TEST_CASE("generate_bam_output_paths: multi input + file path falls back to parent dir",
          "[output_writer]") {
    // Regression: previously produced "out.bam/a_..." for multiple inputs.
    REQUIRE(OutputWriter::generate_bam_output_paths(
        {"a.bam", "b.bam"}, "results/out.bam", "_EASTR_filtered") ==
        std::vector<std::string>{"results/a_EASTR_filtered.bam",
                                 "results/b_EASTR_filtered.bam"});
}
