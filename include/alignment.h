/**
 * Ravi Gaddipati
 * Jan 27, 2016
 * rgaddip1@jhu.edu
 *
 * Provides tools to interact with Alignments.
 *
 * alignment.h
 */

#ifndef VARGAS_ALIGNMENT_H
#define VARGAS_ALIGNMENT_H

#include <vector>
#include "readsource.h"
#include "utils.h"
#include "simdpp/simd.h"
#include "readfile.h"
#include "graph.h"


namespace vargas {

/**
 * Stores a Read and associated alignment information. Positions are
 * 1 indexed.
 * @param read Read that was aligned
 * @param opt_score Best alignment score
 * @param opt_align_end Best alignment position, indexed w.r.t the last base of read
 * @param opt_count Number of alignments that tied for the best score
 * @param sub_score Second best alignment score
 * @param sub_align_end Second best alignment position
 * @param sub_count Number of second-best alignments
 * @param corflag 0 if alignment matches read origin, 1 of second best, 2 otherwise
 */
  struct Alignment {
      Read read;

      // Optimal score information
      uint16_t opt_score;
      int32_t opt_align_end;
      int32_t opt_count;

      // Second best score information
      uint16_t sub_score;
      int32_t sub_align_end;
      int32_t sub_count;

      int8_t corflag;


      Alignment() : opt_score(0), opt_align_end(-1), opt_count(-1), sub_score(0),
                    sub_align_end(-1), sub_count(-1), corflag(-1) { }

      /**
       * Create an alignment with a read and associated meta information.
       * If meta information does not match the expected format, return after
       * populating the raw read.
       * @param line Read
       */
      Alignment(std::string line) {
          Alignment();
          std::vector<std::string> splitLine = split(line, '#');
          // We have the read string
          if (splitLine.size() > 0) {
              this->read.read = splitLine[0];
          }

          // We have meta info
          if (splitLine.size() > 1) {
              splitLine = split(splitLine[1], ',');
              if (splitLine.size() != 12) {
                  // Unexpected format
                  return;
              }

              this->read.end_pos = std::stoi(splitLine[0]);
              this->read.indiv = std::stoi(splitLine[1]);
              this->read.sub_err = std::stoi(splitLine[2]);
              this->read.var_nodes = std::stoi(splitLine[3]);
              this->read.var_bases = std::stoi(splitLine[4]);

              this->opt_score = (uint16_t) std::stoi(splitLine[5]);
              this->opt_align_end = std::stoi(splitLine[6]);
              this->opt_count = std::stoi(splitLine[7]);
              this->sub_score = (uint16_t) std::stoi(splitLine[8]);
              this->sub_align_end = std::stoi(splitLine[9]);
              this->sub_count = std::stoi(splitLine[10]);
              this->corflag = (int8_t) std::stoi(splitLine[11]);
          }
      }

  };

/**
 * Print the alignment to os. This ordering matches the way the alignment is parsed
 * from a string.
 * read_str,opt_score, opt_alignment_end, opt_cout, sub_score, sub_alignment_end, sub_count, corflag
 * Where corflag indicates if the best alignment matched the read origin position (if known) [0], subopt
 * score [1] or other [2].
 * @param os Output stream
 * @param an Alignment output
 */
  inline std::ostream &operator<<(std::ostream &os, const Alignment &a) {
      os << a.read
          << ',' << a.opt_score
          << ',' << a.opt_align_end
          << ',' << a.opt_count
          << ',' << a.sub_score
          << ',' << a.sub_align_end
          << ',' << a.sub_count
          << ',' << int32_t(a.corflag);
      return os;
  }

  /**
   * Aligns a read batch to a reference sequence.
   * All template arguments must match the ReadBatch type.
   * @param read_len maximum read length. Shorter reads are padded by ReadBatch.
   * @param read_len maximum length of the read
   * @param num_reads max number of reads. If a non-default T is used, this should be set to
   *    SIMDPP_FAST_T_SIZE where T corresponds to the width of T. For ex. Default T=simdpp::uint8 uses
   *    SIMDPP_FAST_INT8_SIZE
   * @param CellType element type. Determines max score range. Default uint8
   */
  template<unsigned int read_len,
      unsigned int num_reads = SIMDPP_FAST_INT8_SIZE,
      template<unsigned int, typename=void> class CellType=simdpp::uint8>
  class Aligner {

      /**
       * Default constructor uses the following score values:
       * Match : 2
       * Mismatch : -2
       * Gap Open : 3
       * Gap Extend : 1
       */
      Aligner() { _init(); }

      /**
       * Set scoring parameters
       * @param match match score
       * @param mismatch mismatch score
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      Aligner(int8_t match,
              int8_t mismatch,
              int8_t open,
              int8_t extend) : _match(match), _mismatch(mismatch), _gap_open(open), _gap_extend(extend) { _init(); }

      /**
       * Set the scoring scheme used for the alignments.
       * @param match match score
       * @param mismatch mismatch score
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      void set_scores(int8_t match,
                      int8_t mismatch,
                      int8_t open,
                      int8_t extend) {
          _match = match;
          _mismatch = mismatch;
          _gap_open = open;
          _gap_extend = extend;
          _set_score_matrix();
      }

      /**
       * Align a batch of reads to a topographical sorting of g.
       * @param reads read batch
       * @param g Graph
       * @return vector of alignments
       */
      std::vector<Alignment> align(const ReadBatch<read_len, num_reads, CellType> &reads,
                                   const Graph &g) {
          return align(reads, g.begin(), g.end());
      }

      /**
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param reads batch of reads
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @return vector of alignments
       */
      std::vector<Alignment> align(const ReadBatch<read_len, num_reads, CellType> &reads,
                                   Graph::FilteringIter begin,
                                   Graph::FilteringIter end) {
          std::vector<Alignment> aligns;
          align(reads, begin, end, aligns);
          return aligns;
      }

      /**
      * Align a batch of reads to a graph range, return a vector of alignments
      * corresponding to the reads.
      * @param reads batch of reads
      * @param begin iterator to beginning of graph
      * @param end iterator to end of graph
      * @param aligns vector of alignments
      */
      void align(const ReadBatch<read_len, num_reads, CellType> &reads,
                 Graph::FilteringIter begin,
                 Graph::FilteringIter end,
                 std::vector<Alignment> &aligns) {

      }

    private:
      int8_t _score_matrix[5][5];
      int8_t _match = 2,
          _mismatch = -2,
          _gap_open = 3,
          _gap_extend = 1;

      void _init() {
          _set_score_matrix();
      }

      //TODO this will break if enum Base changes to non-0 start
      inline void _set_score_matrix() {
          _score_matrix[Base::A][Base::A] = _match;
          _score_matrix[Base::A][Base::C] = _mismatch;
          _score_matrix[Base::A][Base::G] = _mismatch;
          _score_matrix[Base::A][Base::T] = _mismatch;
          _score_matrix[Base::A][Base::N] = 0;

          _score_matrix[Base::C][Base::A] = _mismatch;
          _score_matrix[Base::C][Base::C] = _match;
          _score_matrix[Base::C][Base::G] = _mismatch;
          _score_matrix[Base::C][Base::T] = _mismatch;
          _score_matrix[Base::C][Base::N] = 0;

          _score_matrix[Base::G][Base::A] = _mismatch;
          _score_matrix[Base::G][Base::C] = _mismatch;
          _score_matrix[Base::G][Base::G] = _match;
          _score_matrix[Base::G][Base::T] = _mismatch;
          _score_matrix[Base::G][Base::N] = 0;

          _score_matrix[Base::T][Base::A] = _mismatch;
          _score_matrix[Base::T][Base::C] = _mismatch;
          _score_matrix[Base::T][Base::G] = _mismatch;
          _score_matrix[Base::T][Base::T] = _match;
          _score_matrix[Base::T][Base::N] = 0;

          _score_matrix[Base::N][Base::A] = 0;
          _score_matrix[Base::N][Base::C] = 0;
          _score_matrix[Base::N][Base::G] = 0;
          _score_matrix[Base::N][Base::T] = 0;
          _score_matrix[Base::N][Base::N] = 0;
      }
  };

}

#endif //VARGAS_ALIGNMENT_H
