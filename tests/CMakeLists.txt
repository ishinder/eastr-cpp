# Tests CMakeLists.txt

add_executable(eastr_tests
    test_junction.cpp
    test_bed_parser.cpp
    test_strand_handling.cpp
    test_path_utils.cpp
    test_input_detect.cpp
    test_alignment_io.cpp
)

target_link_libraries(eastr_tests PRIVATE
    eastr_lib
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(eastr_tests)

# CLI smoke tests: the binary wiring (incl. --version exiting 0 without --reference).
add_test(NAME cli_version COMMAND eastr --version)
add_test(NAME cli_help COMMAND eastr --help)
