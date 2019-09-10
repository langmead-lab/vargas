/**
 * Ravi Gaddipati
 * Jan 10, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Program CL parsing and scoring structs.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#ifndef VARGAS_SCORING_H
#define VARGAS_SCORING_H

#include "utils.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <string>
#include <sstream>
#include <stdexcept>

namespace vargas {

  /*
   * @brief
   * Marks a forward strand or a reverse complement strand
   */
  enum class Strand {FWD, REV};

  using rg::pos_t;
  using target_t = std::pair<Strand, pos_t>;

  /**
   * @brief
   * Aligner scoring parameters
   */
  struct ScoreProfile {
      ScoreProfile() = default;

      /**
       * @param match Match bonus
       * @param mismatch Mismatch penalty
       * @param gopen Read and ref gap open penalty
       * @param gext Read and ref gap extension penalty
       */
      ScoreProfile(uint8_t match, uint8_t mismatch, uint8_t gopen, uint8_t gext) :
      match(match), mismatch_min(mismatch), mismatch_max(mismatch),
      read_gopen(gopen), read_gext(gext), ref_gopen(gopen), ref_gext(gext), ambig(0) {}

      /**
       * @param match Match bonus
       * @param mismatch Mismatch penalty
       * @param rd_gopen Read gap open penalty
       * @param rd_gext Read gap extension penalty
       * @param ref_gopen Ref gap open penalty
       * @param ref_gext Ref gap extension penalty
       */
      ScoreProfile(unsigned match, unsigned mismatch,
                   unsigned rd_gopen, unsigned rd_gext,
                   unsigned ref_gopen, unsigned ref_gext) :
      match(match), mismatch_min(mismatch), mismatch_max(mismatch),
      read_gopen(rd_gopen), read_gext(rd_gext), ref_gopen(ref_gopen), ref_gext(ref_gext), ambig(0) {}

      /**
       * @brief
       * Get a mismatch penalty from a quality value
       * @param c phred value
       * @return penalty
       */
      unsigned penalty(char c) const {
          return mismatch_min + (unsigned) std::floor( (mismatch_max-mismatch_min) * (std::min<float>(c, 40)/40));
      }

      unsigned
      match = 2, /**< Match bonus */
      mismatch_min = 2,
      mismatch_max = 2, /**< Mismatch penalty */
      read_gopen = 3, /**< Read gap open penalty */
      read_gext = 1, /**< Read gap extension penalty */
      ref_gopen = 3, /**< Ref gap open penalty */
      ref_gext = 1, /**< Ref gap extension penalty */
      ambig = 0; /**< Ambigious base penalty */

      bool end_to_end = false; /**< End to end alignment */

      std::string to_string() const;

  };

  /**
   * @brief
   * Aligner results
   * 1 based coords.
   */
  struct Results {
      std::vector<pos_t> max_pos, sub_pos, max_last_pos, sub_last_pos, waiting_pos, waiting_last_pos;
      std::vector<unsigned> max_count, sub_count;

      std::vector<int> max_score; /**< Best scores */
      std::vector<int> sub_score; /**< Second best scores */

      std::vector<Strand> max_strand;
      std::vector<Strand> sub_strand;

      ScoreProfile profile;

      size_t size() const {
          return max_pos.size();
      }

      /**
       * @brief
       * Resize all result vectors.
       */
      void resize(size_t size);

  };


  const std::vector<std::string> supported_pgid = {"bowtie2", "bwa", "hisat2"};

  std::vector<std::string> tokenize_cl(std::string cl);

  ScoreProfile bwt2(const std::string &cl);

  ScoreProfile bwa_mem(const std::string &cl);

  ScoreProfile program_profile(const std::string &cl);

}

#endif //VARGAS_SCORING_H
