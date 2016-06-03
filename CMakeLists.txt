cmake_minimum_required(VERSION 2.8)

project("vargas")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(libsimdpp libb64/include VByte/include doctest include)

set(MAIN_SOURCES
        src/getopt_pp.cpp       # Argument parsing
        include/getopt_pp.h
        src/utils.cpp           # Global functions
        include/utils.h
        src/graph.cpp           # Graph representation
        include/graph.h
        src/vcfstream.cpp       # VCF file wrapper
        include/vcfstream.h
        # src/readsim.cpp         # Graph read simulator
        # include/readsim.h
        include/readsource.h    # Abstract read source class
        src/readfile.cpp        # Readfile wrapper
        include/readfile.h
        src/xcoder.cpp          # Compressor and base64 encoder for lists of unsigned ints
        include/xcoder.h
        include/alignment.h
        src/samfile.cpp
        include/samfile.h
        include/fasta.h
        include/varfile.h include/dyn_bitset.h)

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
        libb64/include/cencode.h)

set(SSEPP_SOURCES # Source files using libsimdpp
    ${GSSW_SOURCES}
        )


include(libsimdpp/cmake/SimdppMultiarch.cmake)
simdpp_get_runnable_archs(RUNNABLE_ARCHS)


# Dynamic dispatch
# Creates a version of each file for each SSE instruction set.
foreach(SSEPP_FILE ${SSEPP_SOURCES})
    simdpp_multiarch(GEN_ARCH_FILES ${SSEPP_FILE} ${RUNNABLE_ARCHS})
endforeach()
add_library(simdpp ${GEN_ARCH_FILES})

include_directories(${PROJECT_SOURCE_DIR}/htslib)
link_directories(${PROJECT_SOURCE_DIR}/htslib)

add_library(b64 ${B64_SOURCES})
add_library(vbyte ${VBYTE_SOURCES})
add_library(varg ${MAIN_SOURCES})
add_library(gssw ${GSSW_SOURCES})

add_executable(vargas include/main.h src/main.cpp doctest/doctest/doctest.h)
target_link_libraries(vargas simdpp varg gssw vbyte b64 hts)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Global Options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4 -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")

# Release options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O2")

#Debug Options
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -pedantic -Wall -Wextra -Wno-unused-function -Wno-sign-compare")

# Rel w/ Deb options
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} -O3")