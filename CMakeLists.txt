cmake_minimum_required(VERSION 2.8)

project("vmatch")

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(MAIN_SOURCES
        src/getopt_pp.cpp
        include/getopt_pp.h
        src/utils.cpp
        include/utils.h
        src/graph.cpp
        include/graph.h
        src/vcfstream.cpp
        include/vcfstream.h
        src/readsim.cpp
        include/readsim.h
        include/readsource.h
        src/readfile.cpp
        include/readfile.h)

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
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O3")

#Debug Options
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -pedantic -Wall -Wextra -Wno-unused-function")

# Tests
IF (EXISTS "test" AND IS_DIRECTORY "test")
    set(TEST_SOURCES test/graphTest.cpp test/utilsTest.cpp test/vcfStreamTest.cpp test/readsimTest.cpp)
    add_subdirectory(test/googletest/googletest)
    enable_testing()
    include_directories(${gtest_SOURCE_DIR}/include/ ${gtest_SOURCE_DIR}/)
    add_executable(vmatchTests ${TEST_SOURCES})
    target_link_libraries(vmatchTests gtest gtest_main)
    target_link_libraries(vmatchTests vmatch_lib gssw_lib)
ENDIF ()

