#include <catch2/catch_test_macros.hpp>

#include "eastr/alignment_io.hpp"
#include "test_temp.hpp"

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

using namespace eastr;
using eastr_test::TempDir;

namespace {

// 200 bp reference contig; reads below are exact substrings, so CRAM stores them
// as matches against the reference and cannot decode them without it.
const std::string kContig =
    "ACGTTGCAAGGCTTACCGATGCATGCCTAGGATCCAATTGGCCTTAAGGCTAGCTAGGATCGATCGTAGCTAGCTA"
    "GGCTTAACCGGTTAAGGCCTTAAGCTAGCTAGCTAGGATCCGATCGATTACGGATCCGATGCATGCAAGCTTGGC"
    "ATGCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCG";

void write_reference(const std::string& path) {
    std::ofstream(path) << ">chr1\n" << kContig << "\n";
    REQUIRE(fai_build(path.c_str()) == 0);
}

// Write a CRAM with `n_reads` 50 bp reads (one spliced) encoded against `ref`.
void write_cram(const std::string& cram, const std::string& ref, int n_reads) {
    const std::string header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:" +
                               std::to_string(kContig.size()) + "\n";
    samFile* out = sam_open(cram.c_str(), "wc");
    REQUIRE(out != nullptr);
    REQUIRE(hts_set_fai_filename(out, ref.c_str()) == 0);
    sam_hdr_t* hdr = sam_hdr_parse(header.size(), header.c_str());
    REQUIRE(hdr != nullptr);
    REQUIRE(sam_hdr_write(out, hdr) == 0);

    kstring_t line = KS_INITIALIZE;
    bam1_t* b = bam_init1();
    for (int i = 0; i < n_reads; ++i) {
        const int pos = 10 + i * 5;  // 0-based
        std::string rec;
        if (i == 0) {
            // Spliced read: 25M 50N 25M
            rec = "spliced\t0\tchr1\t" + std::to_string(pos + 1) + "\t60\t25M50N25M\t*\t0\t0\t" +
                  kContig.substr(pos, 25) + kContig.substr(pos + 75, 25) + "\t" +
                  std::string(50, 'I');
        } else {
            rec = "r" + std::to_string(i) + "\t0\tchr1\t" + std::to_string(pos + 1) +
                  "\t60\t50M\t*\t0\t0\t" + kContig.substr(pos, 50) + "\t" + std::string(50, 'I');
        }
        line.l = 0;
        kputs(rec.c_str(), &line);
        REQUIRE(sam_parse1(&line, hdr, b) == 0);
        REQUIRE(sam_write1(out, hdr, b) >= 0);
    }
    bam_destroy1(b);
    ks_free(&line);
    sam_hdr_destroy(hdr);
    REQUIRE(sam_close(out) == 0);
}

// Read every record; returns the count, or -1 if decoding fails partway.
int count_records(samFile* in) {
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) return -1;
    bam1_t* b = bam_init1();
    int n = 0, r;
    while ((r = sam_read1(in, hdr, b)) >= 0) ++n;
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    return r == -1 ? n : -1;  // -1 is EOF; < -1 is a decode error
}

// Keep htslib from finding the reference any other way (REF_PATH/REF_CACHE
// lookup, or downloading it by MD5 from the EBI server).
void disable_reference_lookup(const TempDir& d) {
    setenv("REF_PATH", d.file("no_ref_path").c_str(), 1);
    setenv("REF_CACHE", d.file("no_ref_cache").c_str(), 1);
}

} // namespace

TEST_CASE("open_alignment decodes CRAM with the -r reference when UR: is gone",
          "[alignment_io][cram]") {
    TempDir d;
    disable_reference_lookup(d);

    // Encode against a reference, then delete it: the CRAM header's UR: now
    // points nowhere, as with a CRAM made on another machine.
    const std::string original_ref = d.file("original.fa");
    const std::string cram = d.file("reads.cram");
    write_reference(original_ref);
    write_cram(cram, original_ref, 20);
    std::filesystem::remove(original_ref);
    std::filesystem::remove(original_ref + ".fai");

    // The same sequence at a different path, as passed to EASTR with -r.
    const std::string ref = d.file("reference.fa");
    write_reference(ref);

    SECTION("with the reference set, every record decodes") {
        set_alignment_reference(ref);
        samFile* in = open_alignment(cram);
        REQUIRE(in != nullptr);
        REQUIRE(count_records(in) == 20);
        sam_close(in);
    }

    SECTION("without it, decoding fails (what EASTR <= 2.1.1 did)") {
        set_alignment_reference("");
        samFile* in = open_alignment(cram);
        REQUIRE(in != nullptr);
        REQUIRE(count_records(in) == -1);
        sam_close(in);
    }

    set_alignment_reference("");
}

TEST_CASE("open_alignment leaves BAM input unaffected by the reference",
          "[alignment_io]") {
    TempDir d;
    const std::string bam = d.file("reads.bam");
    const std::string header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
    samFile* out = sam_open(bam.c_str(), "wb");
    REQUIRE(out != nullptr);
    sam_hdr_t* hdr = sam_hdr_parse(header.size(), header.c_str());
    REQUIRE(sam_hdr_write(out, hdr) == 0);
    sam_hdr_destroy(hdr);
    REQUIRE(sam_close(out) == 0);

    set_alignment_reference(d.file("does_not_exist.fa"));
    samFile* in = open_alignment(bam);
    REQUIRE(in != nullptr);
    REQUIRE(count_records(in) == 0);
    sam_close(in);
    set_alignment_reference("");
}
