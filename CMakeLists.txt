project("vmatch")

set(MAIN_SOURCES
        main.h
        main.cpp        # Main entry point, launches below functions
        getopt_pp.cpp   # Parser for program arguments
        getopt_pp.h
        buildMain.cpp   # Generate a build file from a ref and VCF
        buildMain.h
        utils.cpp       # Export graph to DOT, split a string, generate a graph from a build file
        utils.h
        alignMain.cpp   # Aligns a read to a graph
        alignMain.h
        simMain.cpp     # Simulates reads from a graph
        simMain.h
        jobMain.cpp     # Generate a HHPC job script (incomplete)
        jobMain.h
)

set(GSSW_SOURCES
        gssw/src/gssw.h
        gssw/src/gssw.c # Main SSW aligner
)

add_executable(vmatch ${MAIN_SOURCES} ${GSSW_SOURCES})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Global Options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4 -std=c++11")

#Debug Options
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -pedantic -Wall -Wextra -Wno-unused-function")