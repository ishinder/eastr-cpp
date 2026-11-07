#include <catch2/catch_test_macros.hpp>

#include "eastr/input_detect.hpp"

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstring>
#include <unistd.h>

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

TEST_CASE("is_bed_file handles a first line longer than NAME_MAX", "[input_detect]") {
    TempFile junk(std::string(400, 'x') + "\n", ".dat");
    REQUIRE(is_bed_file(junk.str()));  // not a path -> BED parser reports the error
}

TEST_CASE("is_bed_file treats a bgzipped file as a single BED", "[input_detect]") {
    TempFile gz(kBamMagic, ".dat");  // gzip magic prefix
    REQUIRE(is_bed_file(gz.str()));
}

TEST_CASE("resolve_input_files returns a single alignment file as-is", "[input_detect]") {
    TempFile bam(kBamMagic, ".dat");
    auto files = resolve_input_files(bam.str(), ListKind::Bam, false);
    REQUIRE(files == std::vector<std::string>{bam.str()});
}

TEST_CASE("resolve_input_files auto-detects a BAM path list", "[input_detect]") {
    // Galaxy-style names: list detection must not depend on ".bam" extensions.
    TempFile a(kBamMagic, ".dat");
    TempFile b(kBamMagic, ".dat");
    TempFile list(a.str() + "\n\n" + b.str() + "\n", ".txt");
    auto files = resolve_input_files(list.str(), ListKind::Bam, false);
    REQUIRE(files == std::vector<std::string>{a.str(), b.str()});
}

TEST_CASE("resolve_input_files with force_list reads the argument as a list",
          "[input_detect]") {
    TempFile a(kBamMagic, ".dat");
    TempFile list(a.str() + "\n", ".txt");
    auto files = resolve_input_files(list.str(), ListKind::Bam, true);
    REQUIRE(files == std::vector<std::string>{a.str()});
}

namespace {
// Message of the exception resolve_input_files throws for `path` (fails if none).
std::string resolve_error(const std::string& path, ListKind kind = ListKind::Bam) {
    try {
        resolve_input_files(path, kind, false);
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    FAIL("expected an exception for " << path);
    return {};
}

bool contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

// A realistic 150 bp record (> NAME_MAX, so a naive fs::exists() would throw).
const std::string kSamRecord =
    "read1\t0\tchr1\t100\t60\t150M\t*\t0\t0\t" + std::string(150, 'A') + "\t" +
    std::string(150, 'I') + "\tNM:i:0\n";
} // namespace

TEST_CASE("resolve_input_files detects a headerless SAM by content", "[input_detect]") {
    // Content-based: works with a Galaxy-style name, not just ".sam".
    TempFile sam(kSamRecord, ".dat");
    std::string msg = resolve_error(sam.str());
    REQUIRE(contains(msg, "looks like a SAM file without a header"));
    REQUIRE(contains(msg, "samtools view -h"));
    REQUIRE_FALSE(contains(msg, "First line:"));  // data line not echoed
}

TEST_CASE("resolve_input_files uses a .sam extension as a hint", "[input_detect]") {
    TempFile sam("not a sam record\n", ".sam");
    REQUIRE(contains(resolve_error(sam.str()), "looks like a SAM file without a header"));
}

TEST_CASE("resolve_input_files flags a .bam/.cram without a valid signature",
          "[input_detect]") {
    TempFile bam("this is not bgzf\n", ".bam");
    std::string msg = resolve_error(bam.str());
    REQUIRE(contains(msg, "has a .bam extension but is not a valid BAM file"));
    REQUIRE(contains(msg, "no BGZF signature"));

    TempFile cram("garbage\n", ".CRAM");  // extension match is case-insensitive
    REQUIRE(contains(resolve_error(cram.str()), "is not a valid CRAM file"));
}

TEST_CASE("resolve_input_files explains an unrecognized tab-separated file",
          "[input_detect]") {
    TempFile junk("a\tb\tc\n", ".dat");
    std::string msg = resolve_error(junk.str());
    REQUIRE(contains(msg, "is not a BAM/CRAM/SAM file"));
    REQUIRE(contains(msg, "--bam takes an alignment file"));
    REQUIRE_FALSE(contains(msg, "First line:"));
}

TEST_CASE("resolve_input_files names a wrong path on the first line", "[input_detect]") {
    TempFile list("/no/such/sample.bam\n", ".txt");
    try {
        resolve_input_files(list.str(), ListKind::Bam, false);
        FAIL("expected an exception");
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        REQUIRE(msg.find("is not a BAM/CRAM/SAM file") != std::string::npos);
        REQUIRE(msg.find("First line: /no/such/sample.bam") != std::string::npos);
    }
}

TEST_CASE("resolve_input_files reports a missing later entry in a list", "[input_detect]") {
    TempFile a(kBamMagic, ".dat");
    TempFile list(a.str() + "\n/no/such/second.bam\n", ".txt");
    try {
        resolve_input_files(list.str(), ListKind::Bam, false);
        FAIL("expected an exception");
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        REQUIRE(msg.find("File not found: /no/such/second.bam") != std::string::npos);
        REQUIRE(msg.find("listed in " + list.str()) != std::string::npos);
    }
}

TEST_CASE("resolve_input_files rejects --bam_list on an actual alignment file",
          "[input_detect]") {
    TempFile bam(kBamMagic, ".bam");
    try {
        resolve_input_files(bam.str(), ListKind::Bam, true);
        FAIL("expected an exception");
    } catch (const std::runtime_error& e) {
        std::string msg = e.what();
        REQUIRE(contains(msg, "is an alignment file, not a list"));
        REQUIRE(contains(msg, "Drop --bam_list"));
        REQUIRE_FALSE(contains(msg, "not a valid BAM"));
    }
}

namespace {
// A scratch directory removed on scope exit.
class TempDir {
public:
    TempDir() {
        path_ = fs::temp_directory_path() /
                ("eastr_dir_" + std::to_string(::getpid()) + "_" + std::to_string(counter_++));
        fs::create_directories(path_);
    }
    ~TempDir() { fs::remove_all(path_); }
    std::string file(const std::string& name) const { return (path_ / name).string(); }
    std::string str() const { return path_.string(); }

private:
    fs::path path_;
    static int counter_;
};
int TempDir::counter_ = 0;

const char* kSamHeader = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";

// Write a real (header-only) BAM/CRAM/SAM via htslib, so format detection is genuine.
void write_alignment(const std::string& path, const char* mode) {
    samFile* out = sam_open(path.c_str(), mode);
    REQUIRE(out != nullptr);
    sam_hdr_t* hdr = sam_hdr_parse(std::strlen(kSamHeader), kSamHeader);
    REQUIRE(hdr != nullptr);
    REQUIRE(sam_hdr_write(out, hdr) == 0);
    sam_hdr_destroy(hdr);
    REQUIRE(sam_close(out) == 0);
}

// Write bgzip-compressed text (gzip magic, but not an alignment file).
void write_bgzf(const std::string& path, const std::string& text) {
    BGZF* fp = bgzf_open(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    REQUIRE(bgzf_write(fp, text.data(), text.size()) == static_cast<ssize_t>(text.size()));
    REQUIRE(bgzf_close(fp) == 0);
}

void write_text(const std::string& path, const std::string& text) {
    std::ofstream(path, std::ios::binary) << text;
}
} // namespace

TEST_CASE("is_alignment_format uses full htslib detection", "[input_detect]") {
    TempDir d;
    write_alignment(d.file("a.dat"), "wb");
    write_bgzf(d.file("genome.fa.gz"), ">chr1\nACGT\n");
    REQUIRE(is_alignment_format(d.file("a.dat")));
    // Same gzip magic as BAM, so the quick check would accept it; htslib does not.
    REQUIRE(is_alignment_file(d.file("genome.fa.gz")));
    REQUIRE_FALSE(is_alignment_format(d.file("genome.fa.gz")));
}

TEST_CASE("resolve_input_files expands a --bam directory to its alignment files",
          "[input_detect]") {
    TempDir d;
    write_alignment(d.file("s2.bam"), "wb");
    write_alignment(d.file("s1.bam"), "wb");
    write_alignment(d.file("s3.sam"), "w");
    write_alignment(d.file("galaxy_dataset.dat"), "wb");  // selected by content
    write_text(d.file("s1.bam.bai"), std::string("BAI\1", 4));  // index: skipped
    write_bgzf(d.file("s2.bam.csi"), "not an alignment");      // gzip but not BAM: skipped
    write_bgzf(d.file("genome.fa.gz"), ">chr1\nACGT\n");     // skipped
    write_text(d.file("notes.txt"), "hello\n");                 // skipped
    write_alignment(d.file(".hidden.bam"), "wb");               // hidden: skipped
    fs::create_directories(d.file("sub"));
    write_alignment(d.file("sub/nested.bam"), "wb");            // not recursive: skipped

    auto files = resolve_input_files(d.str(), ListKind::Bam, false);
    REQUIRE(files == std::vector<std::string>{d.file("galaxy_dataset.dat"), d.file("s1.bam"),
                                              d.file("s2.bam"), d.file("s3.sam")});
}

TEST_CASE("resolve_input_files rejects unusable directories", "[input_detect]") {
    TempDir d;
    write_text(d.file("notes.txt"), "hello\n");
    REQUIRE(contains(resolve_error(d.str()), "no BAM/CRAM/SAM files"));

    write_alignment(d.file("s1.bam"), "wb");
    // --bam_list with a directory
    try {
        resolve_input_files(d.str(), ListKind::Bam, true);
        FAIL("expected an exception");
    } catch (const std::runtime_error& e) {
        REQUIRE(contains(e.what(), "Drop --bam_list"));
    }
    // --bed with a directory is not supported
    REQUIRE(contains(resolve_error(d.str(), ListKind::Bed), "is a directory. --bed takes"));
}
