#include <catch2/catch_test_macros.hpp>

#include "eastr/bed_parser.hpp"

#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;
using namespace eastr;

class TempBedFile {
public:
    TempBedFile(const std::string& content) {
        path_ = fs::temp_directory_path() / ("test_" + std::to_string(counter_++) + ".bed");
        std::ofstream file(path_);
        file << content;
    }

    ~TempBedFile() {
        fs::remove(path_);
    }

    const fs::path& path() const { return path_; }

private:
    fs::path path_;
    static int counter_;
};

int TempBedFile::counter_ = 0;

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

TEST_CASE("BED parser - skip comments and track", "[bed_parser]") {
    TempBedFile bed(
        "# This is a comment\n"
        "track name=junctions\n"
        "chr1\t100\t200\tJUNC1\t5\t+\n"
    );

    auto junctions = BedParser::parse(bed.path().string());

    REQUIRE(junctions.size() == 1);
}

TEST_CASE("BED parser - empty lines", "[bed_parser]") {
    TempBedFile bed(
        "chr1\t100\t200\tJUNC1\t5\t+\n"
        "\n"
        "chr1\t300\t400\tJUNC2\t10\t-\n"
        "\n"
    );

    auto junctions = BedParser::parse(bed.path().string());

    REQUIRE(junctions.size() == 2);
}
