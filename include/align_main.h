/**
 * Ravi Gaddipati
 * Dec 23, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Main aligner interface
 *
 * @file
 */

#ifndef VARGAS_ALIGN_MAIN_H
#define VARGAS_ALIGN_MAIN_H

// Assign all reads without a read group to a RG:ID of:
#define UNGROUPED_READGROUP "VA-UNGROUPED"

// Output alignment tags
#define ALIGN_SAM_MAX_POS_TAG "mp"
#define ALIGN_SAM_SUB_POS_TAG "sp"
#define ALIGN_SAM_MAX_SCORE_TAG "ms"
#define ALIGN_SAM_SUB_SCORE_TAG "ss"
#define ALIGN_SAM_MAX_COUNT_TAG "mc"
#define ALIGN_SAM_SUB_COUNT_TAG "sc"
#define ALIGN_SAM_COR_FLAG_TAG "cf"

#include "cxxopts.hpp"

/**
 * Align given reads to specified target graphs.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int align_main(int argc, char *argv[]);

void align_help(const cxxopts::Options &opts);


#endif //VARGAS_ALIGN_MAIN_H
