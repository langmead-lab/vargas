/**
 * @author Ravi Gaddipati
 * @date June 26, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Interface for simulating and aligning reads from/to a DAG.
 *
 * @file
 */

#ifndef VARGAS_MAIN_H
#define VARGAS_MAIN_H


#include "alignment.h"
#include "graph.h"
#include "readfile.h"
#include "sim.h"

/**
 * Simulate reads from given graph definitions.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int sim_main(const int argc, const char *argv[]);

/**
 * Align given reads to specified target graphs.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int align_main(const int argc, const char *argv[]);

/**
 * Define graphs from a FASTA and a VCF/BCF file. Allows graphs
 * to remain consistent between simulating and aligning steps.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int define_main(const int argc, const char *argv[]);

/**
 * Profile aligner and graph construction.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int profile(const int argc, const char *argv[]);


/**
 * Aligns the given vector of reads to the given graph,
 * using the provided aligners.
 * @param label Label to prepend to output
 * @param subgraph Graph to align to
 * @param reads Reads to align
 * @param aligners Use the given aligners. Size should be equal to number of threads
 * @param out Stream to output result to
 * @param threads number of execution threads.
 */
void align_to_graph(std::string label,
                    const Vargas::Graph &subgraph,
                    const std::vector<Vargas::Read> &reads,
                    const std::vector<std::shared_ptr<Vargas::Aligner<>>> &aligners,
                    std::ostream &out,
                    unsigned int threads);


// Menus
void main_help();
void profile_help();
void align_help();
void sim_help();
void define_help();
#endif //VARGAS_MAIN_H
