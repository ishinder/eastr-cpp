# Tests CMakeLists.txt

add_executable(eastr_tests
    test_main.cpp
    test_junction.cpp
    test_bed_parser.cpp
    test_strand_handling.cpp
)

target_link_libraries(eastr_tests PRIVATE
    eastr_lib
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(eastr_tests)

# Standalone alignment score test (for comparison with Python)
add_executable(test_alignment_score test_alignment_score.cpp)
target_link_libraries(test_alignment_score PRIVATE eastr_lib)
