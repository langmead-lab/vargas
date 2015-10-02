cmake_minimum_required(VERSION 2.6)
project ("vmatch")

set(MAIN_SOURCES
main.h
main.cpp
getopt_pp.cpp
getopt_pp.h
        buildMain.cpp
        buildMain.h utils.cpp utils.h alignMain.cpp alignMain.h simMain.cpp simMain.h jobMain.cpp jobMain.h)

set(GSSW_SOURCES
gssw/src/gssw.h
gssw/src/gssw.c
)

add_executable(vmatch ${MAIN_SOURCES} ${GSSW_SOURCES})
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

add_definitions("-msse4 -std=c++11")