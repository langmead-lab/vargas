/**
 * @author Ravi Gaddipati
 * @date June 26, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Interface for simulating and aligning reads from/to a DAG.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#ifndef VARGAS_MAIN_H
#define VARGAS_MAIN_H

// Version
#define VARGAS_VERSION "dev-" __DATE__

#include "graph.h"
#include "sim.h"
#include "cxxopts.hpp"

#include <stdexcept>

/**
 * Simulate reads from given graph definitions.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int sim_main(int argc, char *argv[]);

/**
 * Define graphs from a FASTA and a VCF/BCF file. Allows graphs
 * to remain consistent between simulating and aligning steps.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int define_main(int argc, char *argv[]);

/**
 * Profile aligner and graph construction.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int profile(int argc, char *argv[]);

/**
 * Extract fields from a SAM file and export them to a CSV file.
 * @param argc CL arg count
 * @param argv CL args
 */
int convert_main(int argc, char **argv);

/**
 * @brief
 * Query sequence files
 * @param argc CL arg count
 * @param argv CL args
 */
int query_main(int argc, char *argv[]);

// Menus
void main_help();
void profile_help(const cxxopts::Options &opts);
void sim_help(const cxxopts::Options &opts);
void define_help(const cxxopts::Options &opts);
void convert_help(const cxxopts::Options &opts);
void query_help(const cxxopts::Options &opts);

#endif //VARGAS_MAIN_H