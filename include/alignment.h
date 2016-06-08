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

#define SIMDPP_ARCH_X86_SSE4_1

#define DEBUG_PRINT_SW 0
#define DEBUG_PRINT_SW_ELEM 0

#include <vector>
#include "readsource.h"
#include "utils.h"
#include "simdpp/simd.h"
#include "readfile.h"
#include "graph.h"
#include "doctest/doctest.h"

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
   * Note: "score" means something that is added, "penalty" refers to something
   * that is substracted. All scores/penalties are provided as positive ints.
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
    public:
      /**
       * Default constructor uses the following score values:
       * Match : 2
       * Mismatch : -2
       * Gap Open : 3
       * Gap Extend : 1
       */
      Aligner() { }

      /**
       * Set scoring parameters
       * @param match match score
       * @param mismatch mismatch penalty
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      Aligner(uint8_t match,
              uint8_t mismatch,
              uint8_t open,
              uint8_t extend) : _match(match), _mismatch(mismatch), _gap_open(open), _gap_extend(extend) { }

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

      void _test_fill_node(const Graph::Node &n,
                           const ReadBatch<read_len, num_reads, CellType> &reads,
                           CellType<num_reads> &max_score,
                           simdpp::uint32<num_reads> &max_pos) {
          _fill_node(n, reads, max_score, max_pos);
      }

    private:
      uint8_t _match = 2,
          _mismatch = 2,
          _gap_open = 3,
          _gap_extend = 1;


      void _fill_node(const Graph::Node &n,
                      const ReadBatch<read_len, num_reads, CellType> &reads,
                      CellType<num_reads> &max_score,
                      simdpp::uint32<num_reads> &max_pos) {

          using namespace simdpp;

          /**
           * Each vector has an 'a' and a 'b' version. Through each row of the
           * matrix fill, their roles are swapped such that one becomes the previous
           * loops data, and the other is filled in.
           * S and D are padded 1 to provide a left column buffer.
           */
          // Score matrix row
          std::vector<CellType<num_reads>, aligned_allocator<CellType<num_reads>, num_reads>>
              Sa(n.seq().size() + 1),
              Sb(n.seq().size() + 1);
          auto S_prev_ptr = &Sa;
          auto S_curr_ptr = &Sb;

          // S_i is a del
          std::vector<CellType<num_reads>, aligned_allocator<CellType<num_reads>, num_reads>>
              Da(n.seq().size() + 1),
              Db(n.seq().size() + 1);
          auto D_prev_ptr = &Da;
          auto D_curr_ptr = &Db;

          // S_i is a ins
          std::vector<CellType<num_reads>, aligned_allocator<CellType<num_reads>, num_reads>>
              Ia(reads.max_len()),
              Ib(reads.max_len());
          auto I_prev_ptr = &Ia;
          auto I_curr_ptr = &Ib;

          // Used for swapping pointers
          std::vector<CellType<num_reads>, aligned_allocator<CellType<num_reads>, num_reads>> *swp_tmp = NULL;

          CellType<num_reads>
              Ceq, // Match score when read_base == ref_base
              Cneq, // mismatch penalty
              tmp;

          max_score = splat(0);

#if DEBUG_PRINT_SW
          std::cout << "   ";
          for (uint32_t col = 1; col < n.seq().size() + 1; ++col) {
              std::cout << num_to_base(static_cast<Base>(n.seq().at(col - 1))) << " ";
          }
#endif

          // For each row of the matrix
          for (size_t row = 0; row < reads.max_len(); ++row) {
              /** Swap the rows we are filling in. The previous row/col becomes what we fill in.
               * S_prev[n] => S(i-1, n)
               * D_prev[n] => D(i-1, n)
               * I_prev[r] => I(r, j-1)
               * where (i,j) is the current cell
               */
              swp_tmp = S_prev_ptr;
              S_prev_ptr = S_curr_ptr;
              S_curr_ptr = swp_tmp;
              const auto &S_prev = *S_prev_ptr;
              auto &S_curr = *S_curr_ptr;

              swp_tmp = D_prev_ptr;
              D_prev_ptr = D_curr_ptr;
              D_curr_ptr = swp_tmp;
              const auto &D_prev = *D_prev_ptr;
              auto &D_curr = *D_curr_ptr;

              swp_tmp = I_prev_ptr;
              I_prev_ptr = I_curr_ptr;
              I_curr_ptr = swp_tmp;
              const auto &I_prev = *I_prev_ptr;
              auto &I_curr = *I_curr_ptr;

#if DEBUG_PRINT_SW
              std::cout << std::endl <<
                  num_to_base(static_cast<Base>(extract<DEBUG_PRINT_SW_ELEM>(reads.at(row)))) << ": ";
#endif

              // Fill in the row. Start at one due to left buffer column
              size_t seq_size = n.seq().size();
              for (uint32_t col = 1; col < seq_size + 1; ++col) {
                  // D(i,j) = D(i-1,j) - gap_extend
                  D_curr[col] = sub_sat(D_prev[col], _gap_extend);
                  // tmp = S(i-1,j) - ( gap_open + gap_extend)
                  tmp = sub_sat(S_prev[col], _gap_extend + _gap_open);
                  // D(i,j) = max{ D(i-1,j) - gap_extend, S(i-1,j) - ( gap_open + gap_extend) }
                  D_curr[col] = max(D_curr[col], tmp);

                  // I(i,j) = I(i,j-1) - gap_extend
                  I_curr[row] = sub_sat(I_prev[row], _gap_extend); // I: I(i,j-1) - gap_extend
                  // tmp = S(i,j-1) - (gap_open + gap_extend)
                  tmp = sub_sat(S_curr[col - 1], _gap_extend + _gap_open);
                  I_curr[row] = max(I_curr[row], tmp);

                  _cmp_read_ref(reads.at(row), n.seq().at(col - 1), Ceq, Cneq); // col - 1 due to column padding
                  S_curr[col] = add_sat(S_prev[col - 1], Ceq); // Add match scores
                  S_curr[col] = sub_sat(S_curr[col], Cneq); // Subtract mismatch scores
                  //S(i,j) = max{ D(i,j), I(i,j), S(i-1,j-1) + C(s,t) }
                  S_curr[col] = max(D_curr[col], S_curr[col]);
                  S_curr[col] = max(I_curr[row], S_curr[col]);

                  // Check for better scores
                  max_score = max(max_score, S_curr[col]);
#if DEBUG_PRINT_SW
                  std::cout << ((int)extract<DEBUG_PRINT_SW_ELEM>(S_curr[col])) << "," << std::flush;
#endif
              }
          }
      }

      /**
       * Compare a ReadBatch position to a reference sequence base.
       * If either base is Base::N, the score/penalty is set to zero.
       * To apply all scores to grid element S[i],
       * S[i] = add_sat(S, Ceq)
       * S[i] = sub_sat(S, Cneq)
       * @param read vector of i'th base of all reads in the ReadBatch
       * @param b reference sequence Base
       * @param Ceq When the read base matches ref seq base, set to match score
       * @param Cneq When the read base doesn't match ref base, set to mismatch penalty
       */
      __attribute__((always_inline))
      inline void _cmp_read_ref(const CellType<num_reads> &read,
                                Base b,
                                CellType<num_reads> &Ceq,
                                CellType<num_reads> &Cneq) {

          using namespace simdpp;

          Ceq = splat(0);
          Cneq = splat(0);

          if (b != Base::N) {
              // Set all mismatching pairs to _mismatch
              auto tmp = cmp_neq(read, b);
              Cneq = tmp & _mismatch;
              // If the read base is Base::N, set to 0 (Ceq)
              tmp = cmp_eq(read, Base::N);
              Cneq = blend(Ceq, Cneq, tmp);

              // b is not N, so all equal bases are valid
              tmp = cmp_eq(read, b);
              Ceq = tmp & _match;
          }
      }

  };

}

TEST_CASE ("Fill Node") {
    std::vector<vargas::Read> reads;
    reads.push_back(vargas::Read("AGTC"));
    reads.push_back(vargas::Read("CATN"));
    reads.push_back(vargas::Read("AGTC"));
    reads.push_back(vargas::Read("CCCC"));

    vargas::ReadBatch<4> rb(reads);

    vargas::Graph::Node n;
    n.set_seq("AGTCATNAGTCCCC");

    simdpp::uint8<16> max;
    simdpp::uint32<16> maxpos;

    vargas::Aligner<4> a;

    a._test_fill_node(n, rb, max, maxpos);
    int i;
}
#endif //VARGAS_ALIGNMENT_H
