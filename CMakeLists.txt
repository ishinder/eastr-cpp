cmake_minimum_required(VERSION 3.16)
project(eastr VERSION 2.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(EASTR_BUILD_TESTS "Build tests" OFF)
option(EASTR_USE_OPENMP "Use OpenMP for parallelization" OFF)

# Find ZLIB (required by htslib)
find_package(ZLIB REQUIRED)

# Find htslib - try multiple methods
set(HTSLIB_FOUND FALSE)

# Method 1: Try pkg-config
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(HTSLIB QUIET htslib>=1.15)
endif()

# Method 2: Try common installation paths
if(NOT HTSLIB_FOUND)
    # Try Homebrew paths on macOS
    find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
        PATHS
            /usr/local/include
            /opt/homebrew/include
            /usr/include
            $ENV{HOME}/miniconda3/include
            $ENV{HOME}/anaconda3/include
    )
    find_library(HTSLIB_LIBRARY
        NAMES hts
        PATHS
            /usr/local/lib
            /opt/homebrew/lib
            /usr/lib
            $ENV{HOME}/miniconda3/lib
            $ENV{HOME}/anaconda3/lib
    )

    if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
        set(HTSLIB_FOUND TRUE)
        set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR})
        set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})
        message(STATUS "Found htslib: ${HTSLIB_LIBRARY}")
    endif()
endif()

if(NOT HTSLIB_FOUND)
    message(FATAL_ERROR "htslib not found. Please install htslib:\n"
        "  macOS: brew install htslib\n"
        "  Ubuntu: apt-get install libhts-dev\n"
        "  conda: conda install -c bioconda htslib")
endif()

# FetchContent for dependencies
include(FetchContent)

# minimap2 - for sequence alignment (same as mappy in Python)
set(MINIMAP2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/minimap2)
if(EXISTS ${MINIMAP2_DIR})
    # Build minimap2 as a static library
    file(GLOB MINIMAP2_SOURCES
        ${MINIMAP2_DIR}/align.c
        ${MINIMAP2_DIR}/bseq.c
        ${MINIMAP2_DIR}/esterr.c
        ${MINIMAP2_DIR}/format.c
        ${MINIMAP2_DIR}/hit.c
        ${MINIMAP2_DIR}/index.c
        ${MINIMAP2_DIR}/jump.c
        ${MINIMAP2_DIR}/kalloc.c
        ${MINIMAP2_DIR}/ksw2_dispatch.c
        ${MINIMAP2_DIR}/ksw2_extd2_sse.c
        ${MINIMAP2_DIR}/ksw2_exts2_sse.c
        ${MINIMAP2_DIR}/ksw2_extz2_sse.c
        ${MINIMAP2_DIR}/ksw2_ll_sse.c
        ${MINIMAP2_DIR}/kthread.c
        ${MINIMAP2_DIR}/lchain.c
        ${MINIMAP2_DIR}/map.c
        ${MINIMAP2_DIR}/misc.c
        ${MINIMAP2_DIR}/options.c
        ${MINIMAP2_DIR}/pe.c
        ${MINIMAP2_DIR}/sdust.c
        ${MINIMAP2_DIR}/seed.c
        ${MINIMAP2_DIR}/sketch.c
        ${MINIMAP2_DIR}/splitidx.c
    )
    add_library(minimap2 STATIC ${MINIMAP2_SOURCES})
    target_include_directories(minimap2 PUBLIC ${MINIMAP2_DIR})
    target_compile_options(minimap2 PRIVATE -msse4.1 -DHAVE_KALLOC)
    target_link_libraries(minimap2 PUBLIC ZLIB::ZLIB pthread)
    set(MINIMAP2_FOUND TRUE)
else()
    message(WARNING "minimap2 not found in external/minimap2. Run: git clone https://github.com/lh3/minimap2.git external/minimap2")
    set(MINIMAP2_FOUND FALSE)
endif()

# CLI11 - command line parser
FetchContent_Declare(
    cli11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.3.2
)
FetchContent_MakeAvailable(cli11)

# OpenMP (optional)
if(EASTR_USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found, enabling parallel loops")
    endif()
endif()

# Main library
add_library(eastr_lib STATIC
    src/types.cpp
    src/junction.cpp
    src/bed_parser.cpp
    src/gtf_parser.cpp
    src/fasta_index.cpp
    src/junction_extractor.cpp
    src/self_aligner.cpp
    src/bowtie2_runner.cpp
    src/spurious_detector.cpp
    src/bam_filter.cpp
    src/output_writer.cpp
    src/thread_pool.cpp
    src/multi_bam_reader.cpp
    src/junction_accumulator.cpp
    src/streaming_spurious_detector.cpp
    src/parallel_bam_extractor.cpp
)

target_include_directories(eastr_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${HTSLIB_INCLUDE_DIRS}
)

target_link_libraries(eastr_lib PUBLIC
    ${HTSLIB_LIBRARIES}
    ZLIB::ZLIB
    CLI11::CLI11
    pthread
)

if(MINIMAP2_FOUND)
    target_link_libraries(eastr_lib PUBLIC minimap2)
    target_compile_definitions(eastr_lib PUBLIC EASTR_USE_MINIMAP2)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(eastr_lib PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(eastr_lib PUBLIC EASTR_USE_OPENMP)
endif()

# Main executable
add_executable(eastr src/main.cpp)
target_link_libraries(eastr PRIVATE eastr_lib)

# Tests
if(EASTR_BUILD_TESTS)
    FetchContent_Declare(
        catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.4.0
    )
    FetchContent_MakeAvailable(catch2)

    enable_testing()
    add_subdirectory(tests)
endif()

# Install
install(TARGETS eastr RUNTIME DESTINATION bin)
install(DIRECTORY include/eastr DESTINATION include)
