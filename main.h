//
// Created by gaddra on 8/6/15.
//

#ifndef VMATCH_H
#define VMATCH_H

#define debug 0

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include "gssw/src/gssw.h"
#include "getopt_pp.h"

/// <summary>
/// Prints node information.
/// </summary>
/// <param name="node">gssw_node to print</param>
void printNode(gssw_node *node);

/// <summary>
/// Splits the specified string, resets elems and returns with split string.
/// </summary>
/// <param name="s">The string</param>
/// <param name="delim">The delimiter</param>
/// <param name="elems">Vector to store results in. Vector is replaced!</param>
/// <returns>Vector of split string.</returns>
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

/// <summary>
/// Builds a previously built graph without the reference or VCF file.
/// </summary>
/// <param name="buildfile">Buildfile generated from a previously built graph</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <returns>Constructed gssw graph.</returns>
gssw_graph* buildGraph(std::string buildfile, int8_t *nt_table, int8_t *mat);


/// <summary>
/// Generates a graph from the given reference and variant file.
/// </summary>
/// <param name="REF">Reference FASTA</param>
/// <param name="VCF">Variant File, uncompressed</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <returns>Constructed gssw_graph</returns>
gssw_graph* generateGraph(std::string REF, std::string VCF, int8_t *nt_table, int8_t *mat, int32_t minpos, int32_t maxpos, std::string outputFile="");


std::string generateRead(gssw_graph &graph, int readLen, float muterr, float indelerr);
#endif //VMATCH_H
