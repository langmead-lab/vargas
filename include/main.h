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
// Program menus

void printBuildHelp();
void printExportHelp();
void printStatHelp();

void main_help();
void profile_help();
void align_help();
void sim_help();
void define_help();

vargas::GraphBuilder load_gdef(std::string file_name,
                               std::unordered_map<std::string, vargas::Graph::Population> &pset);

void align_to_graph(std::string label,
                    const vargas::Graph &subgraph,
                    const std::vector<vargas::Read> &reads,
                    const std::vector<std::shared_ptr<vargas::Aligner<>>> &aligners,
                    std::ostream &out,
                    int threads);

int profile(const int argc, const char *argv[]);
#endif //VARGAS_MAIN_H
