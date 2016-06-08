cmake_minimum_required(VERSION 2.8)

project("vargas")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(libsimdpp doctest include)

set(OPT_PARSE
        src/getopt_pp.cpp       # Argument parsing
        include/getopt_pp.h
        )

set(FILE_IO
        include/fasta.h
        include/varfile.h
        src/readfile.cpp
        include/readfile.h
        include/readsource.h
        )

set(MAIN_SOURCES
        src/utils.cpp           # Global functions
        include/utils.h
        src/graph.cpp           # Graph representation
        include/graph.h
        include/dyn_bitset.h)

set(SIMDPP_SOURCES # Source files using libsimdpp
        include/readsource.h    # Abstract read source class
        include/alignment.h
        )

include_directories(${PROJECT_SOURCE_DIR}/htslib)
link_directories(${PROJECT_SOURCE_DIR}/htslib)

add_library(varg ${MAIN_SOURCES})
add_library(opt ${OPT_PARSE})
add_library(fio ${FILE_IO})

add_executable(vargas include/main.h src/main.cpp)
target_link_libraries(vargas varg hts opt fio)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pg -pedantic -Wall -Wextra -Wno-unused-function -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -O2")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} -O2")