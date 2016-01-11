cmake_minimum_required(VERSION 2.8)

project("vargas")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(libb64/include VByte/include)

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
        src/xcoder.cpp          # Compressor and base64 encoder for lists of unsigned ints
        include/xcoder.h
        include/loadfile.h
        src/loadfile.cpp)

set(GSSW_SOURCES
        gssw/src/gssw.h
        gssw/src/gssw.c # Main SSW aligner
        )

set(VBYTE_SOURCES # Vbyte compression/decompression
        VByte/src/unit.c
        VByte/src/varintdecode.c
        VByte/include/varintdecode.h
        VByte/src/varintencode.c
        VByte/include/varintencode.h
        )

set(B64_SOURCES # B64 encoder/decoder
        libb64/src/cdecode.c
        libb64/src/cencode.c
        libb64/include/encode.h
        libb64/include/decode.h
        libb64/include/cdecode.h
        libb64/include/cencode.h
        )

add_library(b64 ${B64_SOURCES})
add_library(vbyte ${VBYTE_SOURCES})
add_library(varg ${MAIN_SOURCES})
add_library(gssw ${GSSW_SOURCES})

add_executable(vargas include/main.h src/main.cpp)
target_link_libraries(vargas varg gssw vbyte b64)

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
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} -O3")


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
include_directories(${gtest_SOURCE_DIR}/include/ ${gtest_SOURCE_DIR}/) # gtest and gtest_main
add_executable(vargasTests ${TEST_SOURCES})
target_link_libraries(vargasTests gtest gtest_main varg gssw b64 vbyte)
