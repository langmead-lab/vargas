/**
 * Ravi Gaddipati
 * November 25, 2015
 * rgaddip1@jhu.edu
 *
 * Interface for simulating and aligning reads from/to a DAG.
 * Uses a modified gssw from Erik Garrison.
 *
 * main.h
 */

#ifndef VARGAS_MAIN_H
#define VARGAS_MAIN_H

#include "alignment.h"
#include "graph.h"
#include "readfile.h"
#include "sim.h"

// Operational modes
int sim_main(const int argc, const char *argv[]);
int align_main(const int argc, const char *argv[]);
int define_main(const int argc, const char *argv[]);

// Menus
void main_help();
void profile_help();
void align_help();
void sim_help();
void define_help();

/**
 * Load a graph definition file.
 * @param file_name gdef file name
 * @param pset Map to populate label:Population filter pairs
 * @return GraphBuilder to build specified base graph
 */
vargas::GraphBuilder load_gdef(std::string file_name,
                               std::unordered_map<std::string, vargas::Graph::Population> &pset);

/**
 * Aligns the given vector of reads to the given graph,
 * using the provided aligners.
 * @param subgraph Graph to align to
 * @param reads Reads to align
 * @param aligners Use the given aligners. Size should be equal to number of threads
 * @param threads number of execution threads.
 */
void align_to_graph(std::string label,
                    const vargas::Graph &subgraph,
                    const std::vector<vargas::Read> &reads,
                    const std::vector<std::shared_ptr<vargas::Aligner<>>> &aligners,
                    std::ostream &out,
                    int threads);

/**
 * Profile aligner and graph construction.
 */
int profile(const int argc, const char *argv[]);
#endif //VARGAS_MAIN_H
