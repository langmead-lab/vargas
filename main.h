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
/// Prints graph to a DOT format.
/// </summary>
/// <param name="graph">gssw_graph to print</param>
/// <param name="output">output file name</param>
void graphToDot(gssw_graph *graph, std::string output);


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
gssw_graph *buildGraph(std::string buildfile, int8_t *nt_table, int8_t *mat);


/// <summary>
/// Generates a graph from the given reference and variant file.
/// </summary>
/// <param name="REF">Reference FASTA</param>
/// <param name="VCF">Variant File, uncompressed</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <param name="inGroup">Percent of indivudals used to generate graph, default all.</param>
/// <returns>Constructed gssw_graph</returns>
gssw_graph *generateGraph(std::string REF, std::string VCF, int8_t *nt_table, int8_t *mat,
                          int32_t minpos, int32_t maxpos, int32_t maxNodeLen, std::string outputFile = "",
                          int32_t inGroup = -1);

/// <summary>
/// Generates a random read from the given graph. Edges are taken at random.
/// </summary>
/// <param name="graph">Graph to simulate a read from.</param>
/// <param name="readLen">Length of the read to simulate</param>
/// <param name="muterr">Mutation error rate</param>
/// <param name="indelerr">Indel error rate</param>
/// <returns>Constructed gssw_graph</returns>
std::string generateRead(gssw_graph &graph, int readLen, float muterr, float indelerr);


/// <summary>
/// Generates stats from a list of alignment files.
/// </summary>
/// <param name="alignfile">alignment file or comma delimited list</param>
/// <param name="tol">Count as match if the alignment is within this tolerance</param>
void stat_main(std::string alignfile, int32_t tol);


/// <summary>
/// Aligns a list of reads to the provided graphs.
/// </summary>
/// <param name="graph">Graph to align to.</param>
/// <param name="NVgraph">Second graph to align to, empty graph if not used.</param>
/// <param name="mat">Score matrix, use gssw_create_score_matrix()</param>
/// <param name="nt_table">base table, construct using gssw_create_nt_table()</param>
/// <param name="gap_open">Gap open cost</param>
/// <param name="gap_extension">Gap extension cost</param>
/// <param name="readfile">List of input reads</param>
/// <param name="alignfile">Alignment output</param>
/// <param name="dual">True if aligning to graph and NVgraph</param>
/// <returns>Constructed gssw_graph</returns>
void align_main(gssw_graph *graph, gssw_graph *NVgraph,
                int8_t *mat, int8_t *nt_table,
                uint8_t gap_open, uint8_t gap_extension,
                std::string readfile, std::string alignfile, bool dual);

#endif //VMATCH_H
