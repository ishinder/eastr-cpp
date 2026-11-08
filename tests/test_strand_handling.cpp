#include <catch2/catch_test_macros.hpp>

#include "eastr/parallel_bam_extractor.hpp"
#include "test_temp.hpp"

#include <fstream>
#include <string>

using namespace eastr;

TEST_CASE("splice strand comes from XS, not the alignment flag", "[junction][strand]") {
    eastr_test::TempDir dir;
    const std::string sam = dir.file("reads.sam");
    std::ofstream out(sam);
    out << "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n";
    const std::string sequence(50, 'A');
    const std::string quality(50, 'I');
    auto write_read = [&](const char* name, int flag, const char* xs) {
        out << name << '\t' << flag << "\tchr1\t101\t60\t25M50N25M\t*\t0\t0\t"
            << sequence << '\t' << quality << xs << '\n';
    };
    write_read("plus", 16, "\tXS:A:+");
    write_read("minus", 0, "\tXS:A:-");
    write_read("unknown", 16, "");
    out.close();

    const auto junctions = ParallelBamExtractor::extract_junctions_from_bam(sam, "sample", 1);
    REQUIRE(junctions.size() == 3);
    for (Strand strand : {Strand::Plus, Strand::Minus, Strand::Unknown}) {
        const auto it = junctions.find(JunctionKey{"chr1", 125, 175, strand});
        REQUIRE(it != junctions.end());
        REQUIRE(it->second.total_score == 1);
    }
}
