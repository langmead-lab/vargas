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
#define ALIGN_SAM_TARGET_SCORE "ts"
#define ALIGN_SAM_SCORE_PROFILE "pr"

#include "cxxopts.hpp"
#include "sam.h"

// Forward decl to prevent main.cpp recompilation for alignment.h changes
namespace vargas {
  class AlignerBase;
  struct ScoreProfile;
}

/**
 * Align given reads to specified target graphs.
 * @param argc command line argument count
 * @param argv command line arguments
 */
int align_main(int argc, char *argv[]);

/**
 * @brief
 * Create a list of alignment jobs.
 * @param reads input read SAM stream
 * @param align_targets List of targets : RG:Subgraph
 * @param read_len If 0, use length of first read
 * @param chunk_size Limit task size to N alignments
 * @return List of jobs of the form <subgraph label, [reads]>
 */
std::vector<std::pair<std::string, std::vector<vargas::SAM::Record>>>
create_tasks(vargas::isam &reads, std::string &align_targets, const size_t chunk_size, size_t &read_len, bool &resized);

/**
 * @brief
 * Create a new aligner with given parameters
 * @param prof Score profile
 * @param node_len Maximum node length of graph
 * @param read_len Read length
 * @param use_wide use 16 bit cell elements instead of 8 bit
 * @param end_to_end End to end alignment
 * @return pointer to new aligner
 */
std::unique_ptr<vargas::AlignerBase> make_aligner(const vargas::ScoreProfile &prof, size_t read_len,
                                                  bool use_wide = false, bool bare = false);

void align_help(const cxxopts::Options &opts);


#endif //VARGAS_ALIGN_MAIN_H
