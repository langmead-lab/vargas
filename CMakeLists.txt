cmake_minimum_required(VERSION 2.6)
project ("VMatch")
add_executable(VMatch main.h main.cpp gssw/src/gssw.c gssw/src/gssw.h)
add_definitions("-msse4")
