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

#define DEBUG_PRINT_SW 0
#define DEBUG_PRINT_SW_ELEM 6
#define __INLINE__ __attribute__((always_inline)) inline

#include <vector>
#include <algorithm>
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

      Alignment(const Read &r,
                uint16_t best_score, int32_t best_pos, int32_t best_count,
                uint16_t sec_score, int32_t sec_pos, int32_t sec_count,
                int8_t flag) : read(r),
                               opt_score(best_score), opt_align_end(best_pos), opt_count(best_count),
                               sub_score(sec_score), sub_align_end(sec_pos), sub_count(sec_count),
                               corflag(flag) { }

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
   * that is subtracted. All scores/penalties are provided as positive ints.
   * Most memory is allocated for the class so it can be reused during alignment.
   * The state of the memory between function calls is unknown, and must be reset.
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
      typedef typename std::vector<CellType<num_reads>> VecType;

    public:
      /**
       * Default constructor uses the following score values:
       * Match : 2
       * Mismatch : -2
       * Gap Open : 3
       * Gap Extend : 1
       */
      Aligner(size_t max_node_len)
          :
          max_pos(std::vector<uint32_t>(num_reads)),
          max_count(std::vector<uint32_t>(num_reads)),
          sub_pos(std::vector<uint32_t>(num_reads)),
          sub_count(std::vector<uint32_t>(num_reads)),
          _max_node_len(max_node_len) { _alloc(); }

      /**
       * Set scoring parameters
       * @param match match score
       * @param mismatch mismatch penalty
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      Aligner(size_t max_node_len,
              uint8_t match,
              uint8_t mismatch,
              uint8_t open,
              uint8_t extend) :
          _match(match), _mismatch(mismatch), _gap_open(open), _gap_extend(extend),
          max_pos(std::vector<uint32_t>(num_reads)),
          max_count(std::vector<uint32_t>(num_reads)),
          sub_pos(std::vector<uint32_t>(num_reads)),
          sub_count(std::vector<uint32_t>(num_reads)),
          _max_node_len(max_node_len) { _alloc(); }

      ~Aligner() {
          _dealloc();
      }

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
          using namespace simdpp;

          max_score = ZERO_CT;
          max_pos = std::vector<uint32_t>(num_reads);
          _seed seed;
          std::unordered_map<uint32_t, _seed> seed_map;
          for (auto gi = begin; gi != end; ++gi) {
              _get_seed(gi.incoming(), seed_map, &seed);
              seed_map.emplace((*gi).id(), _fill_node(*gi, reads, &seed));
          }

          aligns.clear();
          // Workaround to avoid loops for the template args
          if (num_reads == 16) {
              aligns.emplace(aligns.end(), reads.get_read(0), extract<0>(max_score), max_pos[0], max_count[0],
                             extract<0>(sub_score), sub_pos[0], sub_count[0], 0);
              aligns.emplace(aligns.end(), reads.get_read(1), extract<1>(max_score), max_pos[1], max_count[1],
                             extract<1>(sub_score), sub_pos[1], sub_count[1], 0);
              aligns.emplace(aligns.end(), reads.get_read(2), extract<2>(max_score), max_pos[2], max_count[2],
                             extract<2>(sub_score), sub_pos[2], sub_count[2], 0);
              aligns.emplace(aligns.end(), reads.get_read(3), extract<3>(max_score), max_pos[3], max_count[3],
                             extract<3>(sub_score), sub_pos[3], sub_count[3], 0);
              aligns.emplace(aligns.end(), reads.get_read(4), extract<4>(max_score), max_pos[4], max_count[4],
                             extract<4>(sub_score), sub_pos[4], sub_count[4], 0);
              aligns.emplace(aligns.end(), reads.get_read(5), extract<5>(max_score), max_pos[5], max_count[5],
                             extract<5>(sub_score), sub_pos[5], sub_count[5], 0);
              aligns.emplace(aligns.end(), reads.get_read(6), extract<6>(max_score), max_pos[6], max_count[6],
                             extract<6>(sub_score), sub_pos[6], sub_count[6], 0);
              aligns.emplace(aligns.end(), reads.get_read(7), extract<7>(max_score), max_pos[7], max_count[7],
                             extract<7>(sub_score), sub_pos[7], sub_count[7], 0);
              aligns.emplace(aligns.end(), reads.get_read(8), extract<8>(max_score), max_pos[8], max_count[8],
                             extract<8>(sub_score), sub_pos[8], sub_count[8], 0);
              aligns.emplace(aligns.end(), reads.get_read(9), extract<9>(max_score), max_pos[9], max_count[9],
                             extract<9>(sub_score), sub_pos[9], sub_count[9], 0);
              aligns.emplace(aligns.end(), reads.get_read(10), extract<10>(max_score), max_pos[10], max_count[10],
                             extract<10>(sub_score), sub_pos[10], sub_count[10], 0);
              aligns.emplace(aligns.end(), reads.get_read(11), extract<11>(max_score), max_pos[11], max_count[11],
                             extract<11>(sub_score), sub_pos[11], sub_count[11], 0);
              aligns.emplace(aligns.end(), reads.get_read(12), extract<12>(max_score), max_pos[12], max_count[12],
                             extract<12>(sub_score), sub_pos[12], sub_count[12], 0);
              aligns.emplace(aligns.end(), reads.get_read(13), extract<13>(max_score), max_pos[13], max_count[13],
                             extract<13>(sub_score), sub_pos[13], sub_count[13], 0);
              aligns.emplace(aligns.end(), reads.get_read(14), extract<14>(max_score), max_pos[14], max_count[14],
                             extract<14>(sub_score), sub_pos[14], sub_count[14], 0);
              aligns.emplace(aligns.end(), reads.get_read(15), extract<15>(max_score), max_pos[15], max_count[15],
                             extract<15>(sub_score), sub_pos[15], sub_count[15], 0);
          }
          else {
              throw std::logic_error("No extraction implemented for num reads per vector = " +
                  std::to_string(num_reads));
          }

          // pop off padded reads
          for (size_t i = num_reads; i > reads.reads().size(); --i) aligns.pop_back();
      }

      /**
       * Only used for testing. returns the max score of the read alignment
       * to the node.
       * @param n Graph::Node to align to
       * @param reads ReadBatch to align
       */
      CellType<num_reads> _test_fill_node(const Graph::Node &n,
                                          const ReadBatch<read_len, num_reads, CellType> &reads,
                                          std::vector<uint32_t> &pos) {
          max_score = ZERO_CT;
          max_pos = std::vector<uint32_t>(num_reads);
          _fill_node(n, reads);
          pos = max_pos;
          return max_score;
      }

    protected:

      /**
       * Ending vectors from a previous node
       * @param S_col last column of score matrix
       * @param I_col last column of I vector
       */
      struct _seed {
          _seed() : S_col(read_len), I_col(read_len) { }
          VecType S_col;
          VecType I_col;
      };

      void _get_seed(const std::vector<uint32_t> &prev_ids,
                     const std::unordered_map<uint32_t, _seed> &seed_map,
                     _seed *seed) {
          using namespace simdpp;

          const _seed *ns;

          try {
              for (size_t i = 0; i < read_len; ++i) {
                  seed->I_col[i] = ZERO_CT;
                  seed->S_col[i] = ZERO_CT;
                  for (uint32_t id : prev_ids) {
                      ns = &seed_map.at(id);
                      seed->I_col[i] = max(seed->I_col[i], ns->I_col[i]);
                      seed->S_col[i] = max(seed->S_col[i], ns->S_col[i]);
                  }
                  }
              }
          catch (std::exception &e) {
              std::cerr << "Invalid node ordering." << std::endl;
          }

      }

      /**
       * deletes allocated matrix filling vectors.
       */
      void _dealloc() {
          if (Sa) delete Sa;
          if (Sb) delete Sb;
          if (Da) delete Da;
          if (Db) delete Db;
          if (Ia) delete Ia;
          if (Ib) delete Ib;

          Sa = nullptr;
          Sb = nullptr;
          Da = nullptr;
          Db = nullptr;
          Ia = nullptr;
          Ib = nullptr;

          S_prev = nullptr;
          S_curr = nullptr;
          D_prev = nullptr;
          D_curr = nullptr;
          I_prev = nullptr;
          I_curr = nullptr;
      }

      /**
       * Allocate S and D vectors. I is determined by template parameter.
       * @param seq length of S and D vectors, i.e. node length.
       */
      void _alloc() {
          Sa = new VecType(_max_node_len);
          Sb = new VecType(_max_node_len);
          Da = new VecType(_max_node_len);
          Db = new VecType(_max_node_len);
          Ia = new VecType(read_len);
          Ib = new VecType(read_len);

          S_prev = Sa->data();
          S_curr = Sb->data();
          D_prev = Da->data();
          D_curr = Db->data();
          I_prev = Ia->data();
          I_curr = Ib->data();
      }

      /**
       * Computes local alignment of the node, with no previous seed.
       * @param n Node to align to
       * @param reads ReadBatch to align
       */
      _seed _fill_node(const Graph::Node &n,
                       const ReadBatch<read_len, num_reads, CellType> &reads) {
          _seed s;
          s = _fill_node(n, reads, &s);
          return s;
      }

      /**
       * Computes local alignment to the node.
       * @param n Node to align to
       * @param reads ReadBatch to align
       * @param seeds seeds from previous nodes
       */
      _seed _fill_node(const Graph::Node &n,
                       const ReadBatch<read_len, num_reads, CellType> &reads,
                       const _seed *s) {

          _seed nxt; // Seed for next node
          auto node_seq = n.seq().data(); // Raw pointer faster than at()
          const CellType<num_reads> *read_ptr = reads.data();
          size_t seq_size = n.seq().size();
          uint32_t node_origin = n.end() - seq_size;

          // top left corner
          _fill_cell_rzcz(read_ptr[0], node_seq[0], s);
          _fill_cell_finish(0, 0, node_origin);

          // top row
          for (uint32_t c = 1; c < seq_size; ++c) {
              _fill_cell_rz(read_ptr[0], node_seq[c], c);
              _fill_cell_finish(0, c, node_origin);
          }

          nxt.S_col[0] = S_curr[seq_size - 1];

          // Rest of the rows
          for (uint32_t r = 1; r < read_len; ++r) {
              // Swap the rows we are filling in. The previous row/col becomes what we fill in.
              swp_tmp = S_prev;
              S_prev = S_curr;
              S_curr = swp_tmp;

              swp_tmp = D_prev;
              D_prev = D_curr;
              D_curr = swp_tmp;

              swp_tmp = I_prev;
              I_prev = I_curr;
              I_curr = swp_tmp;

              // first col
              _fill_cell_cz(read_ptr[r], node_seq[0], r, s);
              _fill_cell_finish(r, 0, node_origin);

              // Inner grid
              for (uint32_t c = 1; c < seq_size; ++c) {
                  _fill_cell(read_ptr[r], node_seq[c], r, c);
                  _fill_cell_finish(r, c, node_origin);
              }

              nxt.S_col[r] = S_curr[seq_size - 1];

          }

          // origin vector of what is now I_curr
          nxt.I_col = (Ia->data() == I_curr) ? *Ia : *Ib;
          return nxt;
      }

      /**
       * Fills the top left cell.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param seed alignment seed from previous node
       */
      __INLINE__
      void _fill_cell_rzcz(const CellType<num_reads> &read_base,
                           const Base &ref,
                           const _seed *s) {
          _D(0, ZERO_CT, ZERO_CT);
          _I(0, s->S_col[0]);
          _M(0, read_base, ref, ZERO_CT);
      }

      /**
       * Fills cells when row is 0.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param col current column in matrix
       */
      __INLINE__
      void _fill_cell_rz(const CellType<num_reads> &read_base,
                         const Base &ref,
                         const uint32_t &col) {
          _D(col, ZERO_CT, ZERO_CT);
          _I(0, S_curr[col - 1]);
          _M(col, read_base, ref, ZERO_CT);
      }

      /**
       * Fills cells when col is 0.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row current row in matrix
       * @param seed alignment seed from previous node
       */
      __INLINE__
      void _fill_cell_cz(const CellType<num_reads> &read_base,
                         const Base &ref,
                         const uint32_t &row,
                         const _seed *s) {
          _D(0, D_prev[0], S_prev[0]);
          _I(row, s->S_col[row]);
          _M(0, read_base, ref, s->S_col[row - 1]);
      }

      /**
       * Fills the current cell.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row current row in matrix
       * @param col current column in matrix
       */
      __INLINE__
      void _fill_cell(const CellType<num_reads> &read_base,
                      const Base &ref,
                      const uint32_t &row,
                      const uint32_t &col) {
          using namespace simdpp;

          _D(col, D_prev[col], S_prev[col]);
          _I(row, S_curr[col - 1]);
          _M(col, read_base, ref, S_prev[col - 1]);
      }

      /**
       * Score if there is a deletion
       * @param col current column
       * @param Dp Previous D value at current col.
       * @param Sp Previous S value at current col.
       */
      __INLINE__
      void _D(const uint32_t &col,
              const CellType<num_reads> &Dp,
              const CellType<num_reads> &Sp) {
          using namespace simdpp;

          // D(i,j) = D(i-1,j) - gap_extend
          // Dp is D_prev[col], 0 for row=0
          D_curr[col] = sub_sat(Dp, _gap_extend);
          // tmp = S(i-1,j) - ( gap_open + gap_extend)
          // Sp is S_prev[col], 0 for row=0
          tmp = sub_sat(Sp, _gap_extend + _gap_open);
          // D(i,j) = max{ D(i-1,j) - gap_extend, S(i-1,j) - ( gap_open + gap_extend) }
          D_curr[col] = max(D_curr[col], tmp);
      }

      /**
       * Score if there is an insertion
       * @param row current row
       * @param Sc Previous S value (cell to the left)
       */
      __INLINE__
      void _I(const uint32_t &row,
              const CellType<num_reads> &Sc) {
          using namespace simdpp;

          // I(i,j) = I(i,j-1) - gap_extend
          I_curr[row] = sub_sat(I_prev[row], _gap_extend); // I: I(i,j-1) - gap_extend
          // tmp = S(i,j-1) - (gap_open + gap_extend)
          // Sc is S_curr[col - 1], seed->S_col[row] for col=0
          tmp = sub_sat(Sc, _gap_extend + _gap_open);
          I_curr[row] = max(I_curr[row], tmp);
      }

      /**
       * Best score if there is a match/mismatch. Uses S_prev.
       * @param col current column
       * @param read read base vector
       * @param ref reference sequence base
       * @param Sp Previous S val at col-1 (upper left cell)
       */
      __INLINE__
      void _M(uint32_t col,
              const CellType<num_reads> &read,
              const Base &ref,
              const CellType<num_reads> &Sp) {
          using namespace simdpp;

          Ceq = ZERO_CT;
          Cneq = ZERO_CT;

          if (ref != Base::N) {
              // Set all mismatching pairs to _mismatch
              tmp = cmp_neq(read, ref);
              Cneq = tmp & _mismatch;
              // If the read base is Base::N, set to 0 (Ceq)
              tmp = cmp_eq(read, Base::N);
              Cneq = blend(Ceq, Cneq, tmp);

              // b is not N, so all equal bases are valid
              tmp = cmp_eq(read, ref);
              Ceq = tmp & _match;
          }

          // Sp is S_prev[col - 1], 0 for row=0
          // Seed->S_col[row - 1] for col=0
          S_curr[col] = add_sat(Sp, Ceq); // Add match scores
          S_curr[col] = sub_sat(S_curr[col], Cneq); // Subtract mismatch scores
      }


      /**
       * Takes the max of D,I, and M vectors and stores the current best score/position
       * @param row current row
       * @param col current column
       * @param n Current node, used to get absolute alignment position
       */
      __INLINE__
      __attribute__((optimize("unroll-loops")))
      void _fill_cell_finish(const uint32_t &row,
                             const uint32_t &col,
                             const uint32_t &node_origin) {
          using namespace simdpp;

          //S(i,j) = max{ D(i,j), I(i,j), S(i-1,j-1) + C(s,t) }
          S_curr[col] = max(D_curr[col], S_curr[col]);
          S_curr[col] = max(I_curr[row], S_curr[col]);

          // Check for new high scores
          tmp = S_curr[col] > max_score;
          if (reduce_or(tmp)) {
              max_score = max(max_score, S_curr[col]);
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set
                  if (*(tmp_ptr + (i * _bit_width))) {
                      max_pos[i] = node_origin + col + 1;
                      max_count[i] = 0;
                  }
              }
          }

          // Check for equal max score. If we set a new high score this will set the count to 1
          tmp = cmp_eq(S_curr[col], max_score);
          if (reduce_or(tmp)) {
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set
                  if (*(tmp_ptr + (i * _bit_width))) ++max_count[i];
              }
          }

          // new second best score
          tmp = S_curr[col] > sub_score;
          if (reduce_or(tmp)) {
              uint32_t curr = node_origin + col + 1;
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set
                  if (*(tmp_ptr + (i * _bit_width))) {
                      // Check if far enough from current best
                      if (max_pos[i] < curr - READ_LEN) {
                          //TODO add -(max_score/sub_score) term
                          ;
                      }
                  }
              }
          }

          // Repeat sub score
          tmp = cmp_eq(S_curr[col], sub_score);
          if (reduce_or(tmp)) {
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set
                  if (*(tmp_ptr + (i * _bit_width))) ++sub_count[i];
              }
          }
      }


    private:
      // Used to store max positions with varying vector sizes
      const uint8_t _bit_width = sizeof(CellType<num_reads>) / num_reads;
      // Zero vector
      const CellType<num_reads> ZERO_CT = simdpp::splat(0);

      uint8_t
          _match = 2,       // Match score, is added
          _mismatch = 2,    // mismatch penalty, is subtracted
          _gap_open = 3,    // gap open penalty, subtracted
          _gap_extend = 1;  // gap extension penalty, subtracted

      /**
       * Each vector has an 'a' and a 'b' version. Through each row of the
       * matrix fill, their roles are swapped such that one becomes the previous
       * loops data, and the other is filled in.
       * S and D are padded 1 to provide a left column buffer.
       *
       * The allocated memory persists between node fills, growing if needed.
       */
      VecType
          *Sa = nullptr, // Matrix row
          *Sb = nullptr,
          *Da = nullptr, // Deletion vector
          *Db = nullptr,
          *Ia = nullptr, // Insertion vector
          *Ib = nullptr;

      /**
       * S_prev[n] => S(i-1, n)
       * D_prev[n] => D(i-1, n)
       * I_prev[r] => I(r, j-1)
       * where (i,j) is the current cell
       */
      CellType<num_reads>
          *S_prev = nullptr,
          *S_curr = nullptr,
          *D_prev = nullptr,
          *D_curr = nullptr,
          *I_prev = nullptr,
          *I_curr = nullptr,
          *swp_tmp;

      CellType<num_reads>
          Ceq, // Match score when read_base == ref_base
          Cneq, // mismatch penalty
          tmp; // temporary for use within functions
      const uint8_t *tmp_ptr = (uint8_t *) &tmp;

      CellType<num_reads> max_score;
      std::vector<uint32_t> max_pos;
      std::vector<uint32_t> max_count;

      CellType<num_reads> sub_score;
      std::vector<uint32_t> sub_pos;
      std::vector<uint32_t> sub_count;

      size_t _max_node_len;

  };

}

TEST_CASE ("Alignment") {

    srand(12345);

        SUBCASE("Node Fill scores") {
        std::vector<vargas::Read> reads;
        reads.push_back(vargas::Read("CGT")); // Score 6, Match
        reads.push_back(vargas::Read("ATAGCCA")); // Score 10, del
        reads.push_back(vargas::Read("ATAACGCCA")); // Score 12, Ins
        reads.push_back(vargas::Read("TACCCCA")); // Score 10, mismatch

        vargas::Graph::Node n;
        n.set_endpos(13);
        n.set_seq("ACGTNATACGCCA");

        vargas::ReadBatch<16> rb(reads);
        vargas::Aligner<16> a(15);
        std::vector<uint32_t> pos;

        auto maxscore = a._test_fill_node(n, rb, pos);
            CHECK((int) simdpp::extract<0>(maxscore) == 6);
            CHECK((int) simdpp::extract<1>(maxscore) == 10);
            CHECK((int) simdpp::extract<2>(maxscore) == 12);
            CHECK((int) simdpp::extract<3>(maxscore) == 10);

            CHECK(pos[0] == 4);
            CHECK(pos[1] == 13);
            CHECK(pos[2] == 13);
            CHECK(pos[3] == 13);

    }

        SUBCASE("Graph Alignment") {
        vargas::Graph::Node::_newID = 0;
        vargas::Graph g;

        /**   GGG
        *    /   \
        * AAA     TTT
        *    \   /
        *     CCC(ref)
        */

        {
            vargas::Graph::Node n;
            n.set_endpos(3);
            n.set_as_ref();
            std::vector<bool> a = {0, 1, 1};
            n.set_population(a);
            n.set_seq("AAA");
            g.add_node(n);
        }

        {
            vargas::Graph::Node n;
            n.set_endpos(6);
            n.set_as_ref();
            std::vector<bool> a = {0, 0, 1};
            n.set_population(a);
            n.set_af(0.4);
            n.set_seq("CCC");
            g.add_node(n);
        }

        {
            vargas::Graph::Node n;
            n.set_endpos(6);
            n.set_not_ref();
            std::vector<bool> a = {0, 1, 0};
            n.set_population(a);
            n.set_af(0.6);
            n.set_seq("GGG");
            g.add_node(n);
        }

        {
            vargas::Graph::Node n;
            n.set_endpos(10);
            n.set_as_ref();
            std::vector<bool> a = {0, 1, 1};
            n.set_population(a);
            n.set_seq("TTTA");
            n.set_af(0.3);
            g.add_node(n);
        }

        g.add_edge(0, 1);
        g.add_edge(0, 2);
        g.add_edge(1, 3);
        g.add_edge(2, 3);
        g.set_popsize(3);

        std::vector<vargas::Read> reads;
        reads.push_back(vargas::Read("CCTT"));
        reads.push_back(vargas::Read("GGTT"));
        reads.push_back(vargas::Read("AAGG"));
        reads.push_back(vargas::Read("AACC"));
        reads.push_back(vargas::Read("AGGGT"));
        reads.push_back(vargas::Read("GG"));
        reads.push_back(vargas::Read("AAATTTA"));
        reads.push_back(vargas::Read("AAAGCCC"));

        vargas::ReadBatch<8> rb(reads);
        vargas::Aligner<8> a(5);

        std::vector<vargas::Alignment> aligns = a.align(reads, g);
            REQUIRE(aligns.size() == 8);

            CHECK(aligns[0].read.read == "CCTT");
            CHECK(aligns[0].opt_score == 8);
            CHECK(aligns[0].opt_align_end == 8);

            CHECK(aligns[1].read.read == "GGTT");
            CHECK(aligns[1].opt_score == 8);
            CHECK(aligns[1].opt_align_end == 8);

            CHECK(aligns[2].read.read == "AAGG");
            CHECK(aligns[2].opt_score == 8);
            CHECK(aligns[2].opt_align_end == 5);

            CHECK(aligns[3].read.read == "AACC");
            CHECK(aligns[3].opt_score == 8);
            CHECK(aligns[3].opt_align_end == 5);

            CHECK(aligns[4].read.read == "AGGGT");
            CHECK(aligns[4].opt_score == 10);
            CHECK(aligns[4].opt_align_end == 7);

            CHECK(aligns[5].read.read == "GG");
            CHECK(aligns[5].opt_score == 4);
            CHECK(aligns[5].opt_align_end == 5);

            CHECK(aligns[6].read.read == "AAATTTA");
            CHECK(aligns[6].opt_score == 8);
            CHECK(aligns[6].opt_align_end == 10);

            CHECK(aligns[7].read.read == "AAAGCCC");
            CHECK(aligns[7].opt_score == 8);
            CHECK(aligns[7].opt_align_end == 6);
    }
}

void node_fill_profile() {
        vargas::Graph::Node n;
        {
            std::stringstream ss;
            for (size_t i = 0; i < 10000; ++i) {
                ss << rand_base();
            }
            n.set_seq(ss.str());
        }

        std::vector<vargas::Read> reads;
        for (size_t i = 0; i < 16; ++i) {
            std::stringstream rd;
            for (size_t r = 0; r < 50; ++r) rd << rand_base();
            reads.push_back(vargas::Read(rd.str()));
        }

        vargas::ReadBatch<50> rb(reads);
    vargas::Aligner<50> a(10000);

        std::vector<uint32_t> pos;

    std::cout << "\n10 Node Fill (" << SIMDPP_FAST_INT8_SIZE << "x 50rdlen x 10000bp):\n\t";
        clock_t start = std::clock();
    for (int i = 0; i < 10; ++i) a._test_fill_node(n, rb, pos);
        std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s\n" << std::endl;
}
#endif //VARGAS_ALIGNMENT_H
