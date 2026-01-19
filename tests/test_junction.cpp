#include <catch2/catch_test_macros.hpp>

#include "eastr/junction.hpp"

using namespace eastr;

TEST_CASE("JunctionKey comparison", "[junction]") {
    JunctionKey key1{"chr1", 100, 200, Strand::Plus};
    JunctionKey key2{"chr1", 100, 200, Strand::Plus};
    JunctionKey key3{"chr1", 100, 200, Strand::Minus};
    JunctionKey key4{"chr1", 100, 300, Strand::Plus};
    JunctionKey key5{"chr2", 100, 200, Strand::Plus};

    SECTION("Equality") {
        REQUIRE(key1 == key2);
        REQUIRE_FALSE(key1 == key3);
        REQUIRE_FALSE(key1 == key4);
        REQUIRE_FALSE(key1 == key5);
    }

    SECTION("Less than ordering") {
        REQUIRE(key1 < key5);  // chr1 < chr2
        REQUIRE(key1 < key4);  // same chrom, same start, 200 < 300
        REQUIRE(key1 < key3);  // same coords, Plus < Minus (by char value)
    }

    SECTION("Length calculation") {
        REQUIRE(key1.length() == 100);
        REQUIRE(key4.length() == 200);
    }
}

TEST_CASE("JunctionKey hash", "[junction]") {
    JunctionKeyHash hasher;

    JunctionKey key1{"chr1", 100, 200, Strand::Plus};
    JunctionKey key2{"chr1", 100, 200, Strand::Plus};
    JunctionKey key3{"chr1", 100, 200, Strand::Minus};

    SECTION("Equal keys have equal hashes") {
        REQUIRE(hasher(key1) == hasher(key2));
    }

    SECTION("Different keys have different hashes") {
        // Note: This is probabilistic, but highly likely
        REQUIRE(hasher(key1) != hasher(key3));
    }
}

TEST_CASE("JunctionData operations", "[junction]") {
    JunctionData data;

    SECTION("Add sample") {
        data.add_sample("sample1", 5);
        REQUIRE(data.samples.size() == 1);
        REQUIRE(data.total_score == 5);

        data.add_sample("sample2", 3);
        REQUIRE(data.samples.size() == 2);
        REQUIRE(data.total_score == 8);
    }

    SECTION("Merge junction data") {
        data.add_sample("sample1", 5);

        JunctionData other;
        other.add_sample("sample2", 3);
        other.add_sample("sample3", 2);

        data.merge(other);
        REQUIRE(data.samples.size() == 3);
        REQUIRE(data.total_score == 10);
    }
}

TEST_CASE("Region string parsing", "[junction]") {
    std::string chrom;
    int64_t start, end;

    SECTION("Valid region string") {
        REQUIRE(parse_region_string("chr1:100-200", chrom, start, end));
        REQUIRE(chrom == "chr1");
        REQUIRE(start == 100);
        REQUIRE(end == 200);
    }

    SECTION("Region with underscore in chrom name") {
        REQUIRE(parse_region_string("chr1_random:500-1000", chrom, start, end));
        REQUIRE(chrom == "chr1_random");
        REQUIRE(start == 500);
        REQUIRE(end == 1000);
    }

    SECTION("Invalid region string - no colon") {
        REQUIRE_FALSE(parse_region_string("chr1-100-200", chrom, start, end));
    }

    SECTION("Invalid region string - no dash") {
        REQUIRE_FALSE(parse_region_string("chr1:100", chrom, start, end));
    }
}

TEST_CASE("Make region string", "[junction]") {
    std::string region = make_region_string("chr1", 100, 200);
    REQUIRE(region == "chr1:100-200");
}
