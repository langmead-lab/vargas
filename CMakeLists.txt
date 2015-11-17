cmake_minimum_required(VERSION 3.3)

project("vmatch")

set(MAIN_SOURCES
        src/getopt_pp.cpp   # Parser for program arguments
        include/getopt_pp.h
        src/buildMain.cpp   # Generate a build file from a ref and VCF
        include/buildMain.h
        src/utils.cpp       # Export graph to DOT, split a string, generate a graph from a build file
        include/utils.h
        src/alignMain.cpp   # Aligns a read to a graph
        include/alignMain.h
        src/simMain.cpp     # Simulates reads from a graph
        include/simMain.h
        src/graph.cpp include/graph.h)

set(GSSW_SOURCES
        gssw/src/gssw.h
        gssw/src/gssw.c # Main SSW aligner
        )

add_library(vmatch_lib ${MAIN_SOURCES})
add_library(gssw_lib ${GSSW_SOURCES})
add_executable(vmatch include/main.h src/main.cpp)
target_link_libraries(vmatch vmatch_lib gssw_lib)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Global Options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")

#Debug Options
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -pedantic -Wall -Wextra -Wno-unused-function")

# Tests
add_subdirectory(test/googletest/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include/ ${gtest_SOURCE_DIR}/)
add_executable(vmatchTests test/graphTest.cpp)
target_link_libraries(vmatchTests gtest gtest_main)
target_link_libraries(vmatchTests vmatch_lib gssw_lib)

