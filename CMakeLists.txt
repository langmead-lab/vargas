cmake_minimum_required(VERSION 2.8)

project("vargas")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(gssw/include libb64/include VByte/include)

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
        src/xcoder.cpp
        include/xcoder.h
        )

set(GSSW_SOURCES
        gssw/src/gssw.h
        gssw/src/gssw.c # Main SSW aligner
        )

set(VBYTE_SOURCES
        VByte/src/unit.c
        VByte/src/varintdecode.c
        VByte/include/varintdecode.h
        VByte/src/varintencode.c
        VByte/include/varintencode.h
        )

set(B64_SOURCES
        libb64/src/cdecode.c
        libb64/src/cencode.c
        libb64/include/encode.h
        libb64/include/decode.h
        libb64/include/cdecode.h
        libb64/include/cencode.h
        )

add_library(b64lib ${B64_SOURCES})
add_library(vbytelib ${VBYTE_SOURCES})
add_library(vargaslib ${MAIN_SOURCES})
add_library(gsswlib ${GSSW_SOURCES})

add_executable(vargas include/main.h src/main.cpp)
target_link_libraries(vargas vargaslib gsswlib vbytelib b64lib)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Global Options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4 -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")

# Release options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O3")

#Debug Options
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -pedantic -Wall -Wextra -Wno-unused-function -Wno-sign-compare")

# Rel w/ Deb options
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -O2")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} -O2")


# Tests
set(TEST_SOURCES
        test/graphTest.cpp
        test/utilsTest.cpp
        test/vcfStreamTest.cpp
        test/readsimTest.cpp
        test/equivilence.cpp
        test/equivilence.h
        test/xcoderTest.cpp
        )
add_subdirectory(test/googletest/googletest)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include/ ${gtest_SOURCE_DIR}/)
add_executable(vargasTests ${TEST_SOURCES})
target_link_libraries(vargasTests gtest gtest_main vargaslib gsswlib b64lib vbytelib)
