cmake_minimum_required(VERSION 2.6)
project ("VMatch")

set(MAIN_SOURCES
main.h
main.cpp
getopt_pp.cpp
getopt_pp.h
)
set(GSSW_SOURCES
gssw/src/gssw.h
gssw/src/gssw.c
)

add_executable(VMatch ${MAIN_SOURCES} ${GSSW_SOURCES})

add_definitions("-msse4")
