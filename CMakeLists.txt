cmake_minimum_required(VERSION 2.8)

project("vargas")

FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

include_directories(libsimdpp doctest/doctest include htslib)
link_directories(${PROJECT_SOURCE_DIR}/htslib)

set(MAIN_SOURCES
        src/getopt_pp.cpp
        src/main.cpp
        src/gdef.cpp
        src/graph.cpp
        src/utils.cpp
        src/varfile.cpp
        src/sim.cpp)

set(HEADERS
        include/alignment.h
        include/dyn_bitset.h
        include/fasta.h
        include/gdef.h
        include/graph.h
        include/getopt_pp.h
        include/main.h
        include/sam.h
        include/sim.h
        include/utils.h
        include/varfile.h)


include(libsimdpp/cmake/SimdppMultiarch.cmake)
simdpp_get_runnable_archs(RUNNABLE_ARCHS)

add_executable(vargas ${MAIN_SOURCES} ${HEADERS})
target_link_libraries(vargas hts)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pedantic -Wall -Wextra -Wno-unused-function -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
