# Available Macros:
# BUILD_PHI If set to 'ON' build for Xeon Phi. Requires Intel compiler.
# CMAKE_BUILD_TYPE One of Release, Debug, RelWithDebInfo

cmake_minimum_required(VERSION 2.8)

project("vargas")

FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(doctest/doctest include htslib cxxopts/src)
link_directories(${PROJECT_SOURCE_DIR}/htslib)

set(MAIN_SOURCES
        src/main.cpp
        src/graph.cpp
        src/utils.cpp
        src/varfile.cpp
        src/sim.cpp
        src/fasta.cpp
        src/sam.cpp
        src/align_main.cpp
        src/scoring.cpp
        src/graphman.cpp)

set(HEADERS
        include/alignment.h
        include/dyn_bitset.h
        include/fasta.h
        include/graph.h
        include/main.h
        include/sam.h
        include/sim.h
        include/utils.h
        include/varfile.h
        include/align_main.h
        include/scoring.h
        include/simd.h)

option(BUILD_PHI "Build for Xeon Phi, requires icpc" OFF)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-inline-max-size -wd597 -qoverride-limits")
endif()
option(BUILD_AVX512 "Build for AVX512BW, requires icpc" OFF)

if(BUILD_PHI)
    message("   Building for Xeon Phi")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xMIC-AVX512 -DVA_SIMD_USE_AVX2")
elseif(BUILD_AVX512)
    message("   Building for AVX512")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512bw -mavx512f -DVA_SIMD_USE_AVX512")
else()
    message("   Building for SSE4.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -DVA_SIMD_USE_SSE")
endif()

add_executable(vargas ${MAIN_SOURCES})
target_link_libraries(vargas hts)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pedantic -Wall -Wextra -Wno-unused-function")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_DEBUG}")
