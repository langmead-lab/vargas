//
// Created by gaddra on 9/12/15.
//

#ifndef VMATCH_BUILDMAIN_H
#define VMATCH_BUILDMAIN_H

#include "getopt_pp.h"
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include "utils.h"
#include <algorithm>

#define debug 0

void build_main(int argc, char *argv[]);

/// <summary>
/// Generates a graph from the given reference and variant file.
/// </summary>
/// <param name="REF">Reference FASTA</param>
/// <param name="VCF">Variant File, uncompressed</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <param name="inGroup">Percent of indivudals used to generate graph, default all.</param>
/// <param name="genComplement">Build the complement graph of buildfile</param>
/// <param name="buildfile">Buildfile of the graph to build a complement of.</param>
/// <returns>Constructed gssw_graph</returns>
void generateGraph(
    std::string REF, std::string VCF,
    int32_t minpos, int32_t maxpos, // Region to parse
    int32_t maxNodeLen,
    int32_t inGroup,
    bool genComplement,
    std::string buildfile);


void printBuildHelp();

#endif //VMATCH_BUILDMAIN_H
