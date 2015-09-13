//
// Created by gaddra on 9/12/15.
//

#ifndef VMATCH_ALIGNMAIN_H
#define VMATCH_ALIGNMAIN_H

#include "getopt_pp.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>
#include "utils.h"

void align_main(int argc, char *argv[]);


/// <summary>
/// Aligns a list of reads to the provided graph.
/// </summary>
/// <param name="graph">Graph to align to.</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="gap_open">Gap open cost</param>
/// <param name="gap_extension">Gap extension cost</param>
void align(gssw_graph *graph,
           int8_t *mat, int8_t *nt_table,
           uint8_t gap_open, uint8_t gap_extension);
void printAlignHelp();

#endif //VMATCH_ALIGNMAIN_H
