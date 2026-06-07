#include <catch2/catch_test_macros.hpp>

#include "eastr/input_detect.hpp"

#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;
using namespace eastr;

namespace {

// Temp file that can hold arbitrary (binary) content, with any extension.
class TempFile {
public:
    TempFile(const std::string& content, const std::string& ext = ".dat") {
        path_ = fs::temp_directory_path() /
                ("eastr_detect_" + std::to_string(counter_++) + ext);
        std::ofstream f(path_, std::ios::binary);
        f.write(content.data(), static_cast<std::streamsize>(content.size()));
    }
    ~TempFile() { fs::remove(path_); }
    std::string str() const { return path_.string(); }

private:
    fs::path path_;
    static int counter_;
};
int TempFile::counter_ = 0;

// Minimal BGZF/gzip header followed by the BAM magic.
const std::string kBamMagic = std::string("\x1f\x8b\x08\x04", 4) + "BAM\1";

} // namespace

TEST_CASE("is_alignment_file detects BAM by magic bytes regardless of extension",
          "[input_detect]") {
    // Galaxy names every dataset "*.dat" (issue #2): detection must not use the name.
    TempFile bam(kBamMagic, ".dat");
    REQUIRE(is_alignment_file(bam.str()));
}

TEST_CASE("is_alignment_file detects CRAM and SAM", "[input_detect]") {
    TempFile cram("CRAM\x03", ".dat");
    REQUIRE(is_alignment_file(cram.str()));

    TempFile sam("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", ".dat");
    REQUIRE(is_alignment_file(sam.str()));
}

TEST_CASE("is_alignment_file treats a path list as not-a-single-file", "[input_detect]") {
    TempFile bam(kBamMagic, ".bam");
    TempFile list(bam.str() + "\n", ".txt");
    REQUIRE_FALSE(is_alignment_file(list.str()));
}

TEST_CASE("is_bed_file detects a BED record regardless of extension", "[input_detect]") {
    TempFile bed("chr1\t100\t200\tJUNC1\t5\t+\n", ".dat");
    REQUIRE(is_bed_file(bed.str()));
    // A real BED record is text, not an alignment.
    REQUIRE_FALSE(is_alignment_file(bed.str()));
}

TEST_CASE("is_bed_file skips comment and track lines", "[input_detect]") {
    TempFile bed("# comment\ntrack name=foo\nchr1\t10\t20\tj\t1\t-\n", ".dat");
    REQUIRE(is_bed_file(bed.str()));
}

TEST_CASE("is_bed_file detects a list of existing paths", "[input_detect]") {
    TempFile bed("chr1\t100\t200\tj\t1\t+\n", ".bed");
    TempFile list(bed.str() + "\n", ".txt");
    REQUIRE_FALSE(is_bed_file(list.str()));
}

TEST_CASE("is_bed_file treats a bgzipped file as a single BED", "[input_detect]") {
    TempFile gz(kBamMagic, ".dat");  // gzip magic prefix
    REQUIRE(is_bed_file(gz.str()));
}
