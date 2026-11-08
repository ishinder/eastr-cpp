#include <catch2/catch_test_macros.hpp>

#include "eastr/bed_parser.hpp"
#include "test_temp.hpp"

#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;
using namespace eastr;

class TempBedFile {
public:
    TempBedFile(const std::string& content) : path_(dir_.file("input.bed")) {
        std::ofstream(path_) << content;
    }

    const fs::path& path() const { return path_; }

private:
    eastr_test::TempDir dir_;
    fs::path path_;
};

TEST_CASE("BED parser - basic parsing", "[bed_parser]") {
    TempBedFile bed(
        "chr1\t100\t200\tJUNC1\t5\t+\n"
        "chr1\t300\t400\tJUNC2\t10\t-\n"
        "chr2\t500\t600\tJUNC3\t3\t.\n"
    );

    auto junctions = BedParser::parse(bed.path().string());

    REQUIRE(junctions.size() == 3);

    JunctionKey key1{"chr1", 100, 200, Strand::Plus};
    REQUIRE(junctions.count(key1) == 1);
    REQUIRE(junctions[key1].first == "JUNC1");
    REQUIRE(junctions[key1].second == 5);

    JunctionKey key2{"chr1", 300, 400, Strand::Minus};
    REQUIRE(junctions.count(key2) == 1);
    REQUIRE(junctions[key2].second == 10);
}

TEST_CASE("BED parser skips comments, tracks and blank lines", "[bed_parser]") {
    TempBedFile bed(
        "# This is a comment\n"
        "track name=junctions\n"
        "\n"
        "chr1\t100\t200\tJUNC1\t5\t+\n"
        "\n"
    );

    auto junctions = BedParser::parse(bed.path().string());

    REQUIRE(junctions.size() == 1);
}
