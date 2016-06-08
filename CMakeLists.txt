cmake_minimum_required(VERSION 2.8)

project("vargas")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(libsimdpp doctest include htslib)
link_directories(${PROJECT_SOURCE_DIR}/htslib)

set(FILE_IO
        src/utils.cpp
        src/readfile.cpp
        src/varfile.cpp
        )

set(MAIN_SOURCES
        src/getopt_pp.cpp
        src/main.cpp
        src/graph.cpp
        )

add_library(varg ${MAIN_SOURCES})
add_library(fio ${FILE_IO})

add_executable(vargas ${MAIN_SOURCES})
target_link_libraries(vargas varg hts fio)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pg -pedantic -Wall -Wextra -Wno-unused-function")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -O2")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} -O2")