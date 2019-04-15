# Available Macros:
# BUILD_PHI If set to 'ON' build for Xeon Phi. Requires Intel compiler.
# CMAKE_BUILD_TYPE One of Release, Debug, RelWithDebInfo

cmake_minimum_required(VERSION 2.8)

project("vargas")

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

option(BUILD_KNL "Use Intel compiler to build for KNL Xeon Phi" OFF)
option(BUILD_SKX "Use Intel compiler to build for SKX Xeon Platinum" OFF)
option(BUILD_AVX2 "Use GCC compiler with AVX2 instructions" OFF)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-inline-max-size -wd597 -qoverride-limits")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -wd597")
endif()

option(BUILD_AVX512BW "Build for AVX512BW" OFF)

if(BUILD_KNL)
    message("   Building for KNL Xeon Phi with AVX2 (Intel)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xMIC-AVX2 -DVA_SIMD_USE_AVX2")
elseif(BUILD_AVX512BW)
    message("   Building for SKX Xeon Platinum with AVX512BW (Intel)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCORE-AVX512 -DVA_SIMD_USE_AVX512")
elseif(BUILD_SKX)
    message("   Building for SKX Xeon Platinum with AVX2 (Intel)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xCore-AVX2 -DVA_SIMD_USE_AVX2")
elseif(BUILD_AVX2)
    message("   Building for AVX2 (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -DVA_SIMD_USE_AVX2")
else()
    message("   Building for SSE4.1 (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -DVA_SIMD_USE_SSE")
endif()

add_executable(vargas ${MAIN_SOURCES})
target_link_libraries(vargas hts)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pedantic -Wall -Wextra -Wno-unused-function")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_DEBUG}")
