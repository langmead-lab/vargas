//
// Created by gaddra on 9/12/15.
//

#ifndef VMATCH_UTILS_H
#define VMATCH_UTILS_H

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include "getopt_pp.h"
#include "gssw/src/gssw.h"

/// <summary>
/// Builds a previously built graph without the reference or VCF file.
/// </summary>
/// <param name="buildfile">Buildfile generated from a previously built graph</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <returns>Constructed gssw graph.</returns>
gssw_graph *buildGraph(std::string buildfile, int8_t *nt_table, int8_t *mat);


/// <summary>
/// Splits the specified string, resets elems and returns with split string.
/// </summary>
/// <param name="s">The string</param>
/// <param name="delim">The delimiter</param>
/// <param name="elems">Vector to store results in. Vector is replaced!</param>
/// <returns>Vector of split string.</returns>
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);


/// <summary>
/// Prints node information.
/// </summary>
/// <param name="node">gssw_node to print</param>
void printNode(gssw_node *node);


/// <summary>
/// Prints graph to a DOT format.
/// </summary>
/// <param name="graph">gssw_graph to print</param>
void graphToDot(gssw_graph *graph);

void exportDOT(int argc, char *argv[]);


#endif //VMATCH_UTILS_H
