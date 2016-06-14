cmake_minimum_required(VERSION 2.8)

project("vargas")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(libsimdpp doctest include htslib)
link_directories(${PROJECT_SOURCE_DIR}/htslib)

set(MAIN_SOURCES
        src/getopt_pp.cpp
        src/main.cpp
        src/graph.cpp
        src/utils.cpp
        src/readfile.cpp
        src/varfile.cpp
        )


include(libsimdpp/cmake/SimdppMultiarch.cmake)
simdpp_get_runnable_archs(RUNNABLE_ARCHS)

# SIMD architecture to use
add_definitions(-DSIMDPP_ARCH_X86_SSE4_1)

add_executable(vargas ${MAIN_SOURCES})
target_link_libraries(vargas hts)

target_compile_definitions(vargas PRIVATE READ_LEN=50)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pg -pedantic -Wall -Wextra -Wno-unused-function")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -O2")
