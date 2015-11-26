cmake_minimum_required(VERSION 2.8)

project("vargas")

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(MAIN_SOURCES
        src/getopt_pp.cpp       # Argument parsing
        include/getopt_pp.h
        src/utils.cpp           # Global functions
        include/utils.h
        src/graph.cpp           # Graph representation
        include/graph.h
        src/vcfstream.cpp       # VCF file wrapper
        include/vcfstream.h
        src/readsim.cpp         # Graph read simulator
        include/readsim.h
        include/readsource.h    # Abstract read source class
        src/readfile.cpp        # Readfile wrapper
        include/readfile.h
        )

set(GSSW_SOURCES
        gssw/src/gssw.h
        gssw/src/gssw.c # Main SSW aligner
        )

add_library(vargaslib ${MAIN_SOURCES})
add_library(gsswlib ${GSSW_SOURCES})
add_executable(vargas include/main.h src/main.cpp)
target_link_libraries(vargas vargaslib gsswlib)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Global Options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")

# Release options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O3")

# Rel w/ Deb options
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} -O2")

#Debug Options
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -pedantic -Wall -Wextra -Wno-unused-function -Wno-sign-compare")

# Tests
set(TEST_SOURCES
        test/graphTest.cpp
        test/utilsTest.cpp
        test/vcfStreamTest.cpp
        test/readsimTest.cpp
        test/equivilence.cpp
        test/equivilence.h)
add_subdirectory(test/googletest/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include/ ${gtest_SOURCE_DIR}/)
add_executable(vargasTests ${TEST_SOURCES})
target_link_libraries(vargasTests gtest gtest_main)
target_link_libraries(vargasTests vargaslib gsswlib)
