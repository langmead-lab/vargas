cmake_minimum_required(VERSION 2.8)

project("vargas")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

include_directories(doctest/doctest include htslib cxxopts/src)
link_directories(${PROJECT_SOURCE_DIR}/htslib)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_MACOSX_RPATH 1)
endif()

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

set(LIB_SOURCES
        src/graph.cpp
        src/utils.cpp
        src/varfile.cpp
        src/sim.cpp
        src/fasta.cpp
        src/sam.cpp
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

option(BUILD_AVX512BW_INTEL "Use Intel compiler to build for AVX512BW" OFF)
option(BUILD_AVX512BW_GCC "Use GCC compiler to build for AVX512BW" OFF)
option(BUILD_AVX2_INTEL "Use Intel compiler to build for AVX2" OFF)
option(BUILD_AVX2_GCC "Use GCC compiler to build for AVX2" OFF)

if(BUILD_AVX512BW_INTEL)
    message("   Building for AVX512BW (Intel)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597 -xCORE-AVX512 -DVA_SIMD_USE_AVX512")
elseif(BUILD_AVX512BW_GCC)
    message("   Building for AVX512BW (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512bw -DVA_SIMD_USE_AVX512")
elseif(BUILD_AVX2_INTEL)
    message("   Building for AVX2 (Intel)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd597 -xCORE-AVX2 -DVA_SIMD_USE_AVX2")
elseif(BUILD_AVX2_GCC)
    message("   Building for AVX2 (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -DVA_SIMD_USE_AVX2")
else()
    message("   Building for SSE4.1 (GCC)")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -DVA_SIMD_USE_SSE")
endif()

add_executable(vargas ${MAIN_SOURCES})
add_library(ibvargas SHARED ${LIB_SOURCES})
add_library(static-vargas STATIC ${LIB_SOURCES})
target_link_libraries(vargas hts)
target_link_libraries(ibvargas hts)
target_link_libraries(static-vargas hts)

IF(CMAKE_BUILD_TYPE MATCHES RELEASE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOCTEST_CONFIG_DISABLE")
ENDIF(CMAKE_BUILD_TYPE MATCHES RELEASE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pedantic -Wall -Wextra -Wno-unused-function -fno-inline")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_DEBUG}")
