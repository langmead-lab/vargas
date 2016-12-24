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

include_directories(libsimdpp doctest/doctest include htslib cxxopts/src)
link_directories(${PROJECT_SOURCE_DIR}/htslib)

set(MAIN_SOURCES
        src/main.cpp
        src/gdef.cpp
        src/graph.cpp
        src/utils.cpp
        src/varfile.cpp
        src/sim.cpp
        src/fasta.cpp
        src/sam.cpp
        src/align_main.cpp)

set(HEADERS
        include/alignment.h
        include/dyn_bitset.h
        include/fasta.h
        include/gdef.h
        include/graph.h
        include/main.h
        include/sam.h
        include/sim.h
        include/utils.h
        include/varfile.h
        include/align_main.h)

option(BUILD_PHI "Build for Xeon Phi, requires icpc" OFF)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -no-inline-max-size -qoverride-limits -qopt-report=4")
endif ()

if (BUILD_PHI)
    message("Building for Xeon Phi")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xMIC-AVX512 -DSIMDPP_ARCH_X86_AVX512F")
else ()
    message("Building for SSE4.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -DSIMDPP_ARCH_X86_SSE4_1")
endif ()

add_executable(vargas ${MAIN_SOURCES})
target_link_libraries(vargas hts)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pedantic -Wall -Wextra -Wno-unused-function")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_DEBUG}")
