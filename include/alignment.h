/**
 * @author Ravi Gaddipati
 * @date June 21, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Provides tools to interact with Alignments.
 *
 * @details
 * Reads are aligned using
 * SIMD vectorized smith-waterman. Reads are grouped into batches, with
 * size determined by the optimal hardware support. The default template
 * parameters support an 8 bit score.
 *
 * @file
 */

#ifndef VARGAS_ALIGNMENT_H
#define VARGAS_ALIGNMENT_H


#define DEBUG_PRINT_SW 1 // Print the full SW matrix for each node aligned
#define DEBUG_PRINT_SW_NUM 0 // Print the matrix for this read number in the alignment group

#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include "readsource.h"
#include "utils.h"
#include "simdpp/simd.h"
#include "readfile.h"
#include "graph.h"
#include "doctest.h"

namespace Vargas {

/**
 * @brief
 * Read alignment storage struct.
 * @details
 * Stores a Read and associated alignment information. Best alignment information as well
 * as the second best alignment information is stored. The alignment closest to the target is kept.
 */
  struct Alignment {
      Read read;

      // Optimal score information
      uint16_t opt_score;
      /**< Read that was aligned. */
      int32_t opt_align_end;
      /**< Best alignment position, indexed w.r.t the last base of read.*/
      int32_t opt_count; /**< Number of alignments that tied for the best score.*/

      // Second best score information
      uint16_t sub_score;
      /**< Second best alignment score.*/
      int32_t sub_align_end;
      /**< Second best alignment position, indexed w.r.t the last base of the read.*/
      int32_t sub_count;
      /**< Number of second-best alignments.*/

      int8_t corflag;
      /**< 1 if Read origin matched best score, 2 for sub score, 0 otherwise.*/

      Alignment() : opt_score(0), opt_align_end(-1), opt_count(-1), sub_score(0),
                    sub_align_end(-1), sub_count(-1) { }

      Alignment(const Read &r,
                uint16_t best_score, int32_t best_pos, int32_t best_count,
                uint16_t sec_score, int32_t sec_pos, int32_t sec_count, uint8_t cor) :
          read(r),
          opt_score(best_score), opt_align_end(best_pos), opt_count(best_count),
          sub_score(sec_score), sub_align_end(sec_pos), sub_count(sec_count), corflag(cor) { }

  };

  /**
   * @brief Main SIMD SW Aligner.
   * @details
   * Aligns a read batch to a reference sequence.
   * Note: "score" means something that is added, "penalty" refers to something
   * that is subtracted. All scores/penalties are provided as positive integers.
   * Most memory is allocated for the class so it can be reused during alignment. To reduce memory usage,
   * the maximum node size can be reduced.
   * @tparam num_reads max number of reads. If a non-default T is used, this should be set to
   *    SIMDPP_FAST_T_SIZE where T corresponds to the width of T. For ex. Default T=simdpp::uint8 uses
   *    SIMDPP_FAST_INT8_SIZE
   * @tparam CellType element type. Determines max score range. Default simdpp::uint8.
   * @tparam NativeT Native version of CellType. Default uint8_t.
   */
  template<unsigned int num_reads = SIMDPP_FAST_INT8_SIZE,
      template<unsigned int, typename=void> class CellType=simdpp::uint8,
      typename NativeT=uint8_t>
  class Aligner {
      typedef typename std::vector<CellType<num_reads>> VecType;

    public:
      /**
       * @brief
       * Default constructor uses the following score values: \n
       * Match : 2 \n
       * Mismatch : -2 \n
       * Gap Open : 3 \n
       * Gap Extend : 1 \n
       * @param max_node_len maximum node length
       * @param len maximum read length
       */
      Aligner(size_t max_node_len, int len)
          :
          read_len(len),
          max_pos(std::vector<uint32_t>(num_reads)),
          max_count(std::vector<uint32_t>(num_reads)),
          sub_pos(std::vector<uint32_t>(num_reads)),
          sub_count(std::vector<uint32_t>(num_reads)),
          _max_node_len(max_node_len) { _alloc(); }

      /**
       * @brief
       * Set scoring parameters.
       * @param max_node_len max node length
       * @param len max read length
       * @param match match score
       * @param mismatch mismatch penalty
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      Aligner(size_t max_node_len,
              int len,
              uint8_t match,
              uint8_t mismatch,
              uint8_t open,
              uint8_t extend) :
          read_len(len),
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
 * @brief
 * Container for a packaged batch of reads.
 * @details
 * Reads are interleaved so each SIMD vector
 * contains bases from all reads, respective to the base number. For example AlignmentGroup[0]
 * would contain the first bases of every read. Short reads or missing reads are padded
 * with Base::N.
 * @tparam num_reads max number of reads. If a non-default T is used, this should be set to
 *    SIMDPP_FAST_T_SIZE where T corresponds to the width of T. For ex. Default T=simdpp::uint8 uses
 *    SIMDPP_FAST_INT8_SIZE
 * @tparam T element type
 */
      class AlignmentGroup {
        public:

          /**
           * @param len maximum read length
           */
          AlignmentGroup(int len) : read_len(len) { }

          /**
           * @brief
           * Read length is set to first read size.
           * @param batch package the given vector of reads. Must be nonempty.
           */
          AlignmentGroup(const std::vector<std::vector<Base>> &batch) : _reads(batch) {
              if (batch.size() == 0) throw std::invalid_argument("Vector of reads must be non-empty.");
              read_len = batch[0].size();
              _package_reads();
          }

          /**
           * @param batch package the given vector of reads. Must be nonempty.
           * @param len max read length
           */
          AlignmentGroup(const std::vector<std::vector<Base>> &batch, int len) : read_len(len), _reads(batch) {
              _package_reads();
          }

          /**
           * @brief
           * Read length is set to first read size.
           * @param batch package the given vector of reads. Must be nonempty.
           */
          AlignmentGroup(const std::vector<std::string> &batch) {
              if (batch.size() == 0) throw std::invalid_argument("Vector of reads must be non-empty.");
              _reads.clear();
              for (auto &seq : batch) _reads.push_back(seq_to_num(seq));
              read_len = batch[0].length();
              _package_reads();
          }

          /**
           * @param batch package the given vector of reads. Must be nonempty.
           * @param len max read length
           */
          AlignmentGroup(const std::vector<std::string> &batch, int len) : read_len(len) {
              _reads.clear();
              for (auto &seq : batch) _reads.push_back(seq_to_num(seq));
              _package_reads();
          }

          /**
           * @param batch load the given vector of reads.
           */
          bool load_reads(const std::vector<std::string> &batch) {
              if (batch.size() == 0) return false;
              _reads.clear();
              for (auto &seq : batch) {
                  _reads.push_back(seq_to_num(seq));
              }
              _package_reads();
              return true;
          }

          /**
           * @param batch load the given vector of reads.
           */
          bool load_reads(const std::vector<std::vector<Base>> &batch) {
              if (batch.size() == 0) return false;
              _reads = batch;
              _package_reads();
              return true;
          }

          /**
           * @brief
           * Return the i'th base of every read in a simdpp vector.
           * @param i base index.
           */
          const CellType<num_reads> &at(int i) const {
              // let vector handle out of range errors
              return _packaged_reads.at(i);
          }

          /**
           * @brief
           * Pointer to raw packaged read data.
           * @return CellType<num_reads> pointer
           */
          const CellType<num_reads> *data() const {
              return _packaged_reads.data();
          }

          /**
           * @brief
           * Non const version of at(i).
           * @param i base index
           */
          CellType<num_reads> &operator[](int i) {
              return _packaged_reads.at(i);
          }

          /**
           * @return max read length. Echos template arg
           */
          size_t max_len() const { return read_len; }

          /**
           * @brief
           * Returns optimal number of reads in a batch based on SIMD architecture.
           * @return batch size.
           */
          size_t batch_size() const { return num_reads; }

          /**
           * @return Reads used to build the batch.
           */
          const std::vector<std::vector<Base>> &reads() const { return _reads; }

          /**
           * @brief
           * Get a read, empty read if out of range.
           * @param i index of read
           * @return Read at i
           */
          std::vector<Base> get_read(int i) const {
              return _reads.at(i);
          }

          /**
           * @brief
           * Get the utilization of the batch capacity. In effect how much
           * padding was used.
           * @return fill, between 0 and 1.
           */
          float fill() const {
              float f = 0;
              for (auto &r : _reads) f += r.size();
              return f / (num_reads * read_len);
          }

          typename std::vector<CellType<num_reads>>::const_iterator begin() const { return _packaged_reads.begin(); }
          typename std::vector<CellType<num_reads>>::const_iterator end() const { return _packaged_reads.end(); }

        private:

          int read_len;

          /**
           * _packaged_reads[i] contains all i'th bases.
           * The length of _packaged_reads is the length of the read,
           * where as the length of _packaged_reads[i] is the number
           * of reads.
           */
          std::vector<CellType<num_reads>> _packaged_reads;

          // Unpackaged reads
          std::vector<std::vector<Base>> _reads;

          /**
           * Interleaves reads so all same-index base positions are in one
           * vector. Empty spaces are padded with Base::N.
           */
          inline void _package_reads() {
              _packaged_reads.resize(read_len);
              if (_reads.size() > num_reads) throw std::range_error("Too many reads for batch size.");

              // allocate memory
              uchar **pckg = (uchar **) malloc(read_len * sizeof(uchar *));
              for (int i = 0; i < read_len; ++i) {
                  pckg[i] = (uchar *) malloc(num_reads * sizeof(uchar));
              }

              // Interleave reads
              // For each read (read[i] is in _packaged_reads[0..n][i]
              for (size_t r = 0; r < _reads.size(); ++r) {

                  if (_reads.at(r).size() > read_len) throw std::range_error("Read too long for batch size.");

                  // Put each base in the appropriate vector element
                  for (size_t p = 0; p < _reads[r].size(); ++p) {
                      pckg[p][r] = _reads[r][p];
                  }

                  // Pad the shorter reads
                  for (size_t p = _reads[r].size(); p < read_len; ++p) {
                      pckg[p][r] = Base::N;
                  }
              }

              // Pad underful batches
              for (size_t r = _reads.size(); r < num_reads; ++r) {
                  for (size_t p = 0; p < read_len; ++p) {
                      pckg[p][r] = Base::N;
                  }
              }

              // Load into vectors
              for (int i = 0; i < read_len; ++i) {
                  _packaged_reads[i] = simdpp::load(pckg[i]);
              }

              // Free memory
              for (int i = 0; i < read_len; ++i) {
                  free(pckg[i]);
              }
              free(pckg);

          }

      };

      /**
       * @brief
       * Struct to return the alignment results
       */
      struct Results {
          Results() : best_pos(num_reads), sub_pos(num_reads), best_score(num_reads), sub_score(num_reads) { }
          std::vector<uint32_t> best_pos;
          /**< Best positions */
          std::vector<uint32_t> sub_pos;
          /**< Second best positions */
          std::vector<uint8_t> best_score;
          /**< Best scores */
          std::vector<uint8_t> sub_score; /**< Second best scores */
      };

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
       * @brief
       * Align a batch of reads to the given graph.
       * @param read_group AlignmentGroup to align
       * @param targets correct positions for reads in read_group
       * @param g Graph
       * @return Results packet
       */
      Results align(const AlignmentGroup &read_group,
                    const std::vector<uint32_t> &targets,
                    const Graph &g) {
          return align(read_group, targets, g.begin(), g.end());
      }

      /**
       * @brief
       * Align a batch of reads to the given graph.
       * @param read_group AlignmentGroup to align
       * @param g Graph
       * @return Results packet
       */
      Results align(const AlignmentGroup &read_group,
                    const Graph &g) {
          std::vector<uint32_t> targets(num_reads);
          return align(read_group, targets, g.begin(), g.end());
      }

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group AlignmentGroup to align
       * @param targets correct positions for reads in read_group
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @return Results packet
       */
      Results align(const AlignmentGroup &read_group,
                    const std::vector<uint32_t> &targets,
                    Graph::FilteringIter begin,
                    Graph::FilteringIter end) {
          Results aligns;
          align_into(read_group, targets, begin, end, aligns);
          return aligns;
      }

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group AlignmentGroup to align
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @return Results packet
       */
      Results align(const AlignmentGroup &read_group,
                    Graph::FilteringIter begin,
                    Graph::FilteringIter end) {
          Results aligns;
          std::vector<uint32_t> targets(num_reads);
          align_into(read_group, targets, begin, end, aligns);
          return aligns;
      }

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group AlignmentGroup to align
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @param aligns Results packet to populate
     */
      void align_into(const AlignmentGroup &read_group,
                      Graph::FilteringIter begin,
                      Graph::FilteringIter end,
                      Results &aligns) {
          std::vector<uint32_t> targets(num_reads);
          align_into(read_group, targets, begin, end, aligns);
      }

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group AlignmentGroup to align
       * @param targets correct positions for reads in read_group
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @param aligns Results packet to populate
       */
      void align_into(const AlignmentGroup &read_group,
                      const std::vector<uint32_t> &targets,
                      Graph::FilteringIter begin,
                      Graph::FilteringIter end,
                      Results &aligns) {
          using namespace simdpp;

          // Reset old info
          max_score = ZERO_CT;
          sub_score = ZERO_CT;
          max_pos = std::vector<uint32_t>(num_reads);
          sub_pos = std::vector<uint32_t>(num_reads);
          _seed seed(read_len);

          std::unordered_map<uint32_t, _seed> seed_map;
          for (auto gi = begin; gi != end; ++gi) {
              _get_seed(gi.incoming(), seed_map, &seed);
              if ((*gi).is_pinched()) seed_map.clear();
              seed_map.emplace((*gi).id(), _fill_node(*gi, read_group, targets, &seed));
          }

          for (int i = 0; i < num_reads; ++i) {
              aligns.best_score[i] = extract(i, max_score);
              aligns.sub_score[i] = extract(i, sub_score);
              aligns.best_pos[i] = max_pos[i];
              aligns.sub_pos[i] = sub_pos[i];
          }
      }

      /**
       * @brief
       * Only used for testing. returns the max score of the read alignment
       * to the node.
       * @param n Node to align to
       * @param read_group AlignmentGroup to align
       * @param pos positions of the max score
       */
      CellType<num_reads> _test_fill_node(const Graph::Node &n,
                                          const AlignmentGroup &read_group,
                                          std::vector<uint32_t> &pos) {
          max_score = ZERO_CT;
          sub_score = ZERO_CT;
          std::vector<uint32_t> targets(num_reads);
          max_pos = std::vector<uint32_t>(num_reads);
          _fill_node(n, read_group, targets);
          pos = max_pos;
          return max_score;
      }

    protected:

      /**
       * @brief
       * Ending vectors from a previous node
       */
      struct _seed {
          _seed(int read_len) : S_col(read_len), I_col(read_len) { }
          VecType S_col;
          /**< Last column of score matrix.*/
          VecType I_col;/**< Last column of I vector.*/
      };

      /**
       * @brief
       * Returns the best seed from all previous nodes.
       * @param prev_ids All nodes preceding current node
       * @param seed_map ID->seed map for all previous nodes
       * @param seed best seed to populate
       */
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
              throw std::logic_error("Unable to get seed, invalid node ordering.");
          }

      }

      /**
       * @brief
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
       * @brief
       * Allocate S and D vectors. I is determined by template parameter.
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
       * @brief
       * Computes local alignment of the node, with no previous seed.
       * @param n Node to align to
       * @param read_group AlignmentGroup to align
       * @param targets correct positions for reads in read_group
       */
      _seed _fill_node(const Graph::Node &n,
                       const AlignmentGroup &read_group,
                       const std::vector<uint32_t> &targets) {
          _seed s(read_len);
          s = _fill_node(n, read_group, targets, &s);
          return s;
      }

      /**
       * @brief
       * Computes local alignment to the node.
       * @param n Node to align to
       * @param read_group AlignmentGroup to align
       * @param targets correct positions for reads in read_group
       * @param s seeds from previous nodes
       */
      _seed _fill_node(const Graph::Node &n,
                       const AlignmentGroup &read_group,
                       const std::vector<uint32_t> &targets,
                       const _seed *s) {

          _seed nxt(read_len); // Seed for next node
          auto node_seq = n.seq().data();
          const CellType<num_reads> *read_ptr = read_group.data();
          size_t seq_size = n.seq().size();
          uint32_t node_origin = n.end() - seq_size;

          #if DEBUG_PRINT_SW
          std::cout << std::endl << "-\t";
          for (auto c : n.seq_str()) std::cout << c << '\t';
          std::cout << std::endl;
          #endif

          // top left corner
          _fill_cell_rzcz(read_ptr[0], node_seq[0], s);
          _fill_cell_finish(0, 0, node_origin, targets);

          // top row
          for (uint32_t c = 1; c < seq_size; ++c) {
              _fill_cell_rz(read_ptr[0], node_seq[c], c);
              _fill_cell_finish(0, c, node_origin, targets);
          }

          #if DEBUG_PRINT_SW
          std::cout << num_to_base((Base) simdpp::extract<DEBUG_PRINT_SW_NUM>(read_ptr[0])) << '\t';
          for (size_t i = 0; i < n.seq().size(); ++i)
              std::cout << (int) simdpp::extract<DEBUG_PRINT_SW_NUM>(S_curr[i]) << '\t';
          std::cout << std::endl;
          #endif

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
              _fill_cell_finish(r, 0, node_origin, targets);

              // Inner grid
              for (uint32_t c = 1; c < seq_size; ++c) {
                  _fill_cell(read_ptr[r], node_seq[c], r, c);
                  _fill_cell_finish(r, c, node_origin, targets);
              }

              nxt.S_col[r] = S_curr[seq_size - 1];

              #if DEBUG_PRINT_SW
              std::cout << num_to_base((Base) simdpp::extract<DEBUG_PRINT_SW_NUM>(read_ptr[r])) << '\t';
              for (size_t i = 0; i < n.seq().size(); ++i)
                  std::cout << (int) simdpp::extract<DEBUG_PRINT_SW_NUM>(S_curr[i]) << '\t';
              std::cout << std::endl;
              #endif

          }

          // origin vector of what is now I_curr
          nxt.I_col = (Ia->data() == I_curr) ? *Ia : *Ib;
          return nxt;
      }

      /**
       * @brief
       * Fills the top left cell.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param s alignment seed from previous node
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
       * @brief
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
       * @brief
       * Fills cells when col is 0.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row current row in matrix
       * @param s alignment seed from previous node
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
       * @brief
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
       * @brief
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
       * @brief
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
       * @brief
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
       * @brief
       * Takes the max of D,I, and M vectors and stores the current best score/position
       * Currently does not support non-deafault template args
       * @param row current row
       * @param col current column
       * @param node_origin Current position, used to get absolute alignment position
       * @param targets Used to mark corflag, indicates 'correct' alignments
       */
      __INLINE__
      void _fill_cell_finish(const uint32_t &row,
                             const uint32_t &col,
                             const uint32_t &node_origin,
                             const std::vector<uint32_t> &targets) {
          using namespace simdpp;

          uint32_t curr = node_origin + col + 1; // absolute position

          //S(i,j) = max{ D(i,j), I(i,j), S(i-1,j-1) + C(s,t) }
          S_curr[col] = max(D_curr[col], S_curr[col]);
          S_curr[col] = max(I_curr[row], S_curr[col]);

          // Check for new high scores
          tmp = S_curr[col] > max_score;
          NativeT max_elem;
          if (reduce_or(tmp)) {
              max_score = max(max_score, S_curr[col]);
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set
                  if (extract(i, tmp)) {
                      max_elem = extract(i, max_score);
                      // If old max is larger than old sub_max, and if its far enough away
                      if (max_elem > extract(i, sub_score) && curr < max_pos[i] - read_len) {
                          insert(max_elem, i, sub_score);
                          sub_pos[i] = max_pos[i];
                          sub_count[i] = max_count[i];
                      }
                      max_pos[i] = curr;
                      max_count[i] = 0;
                  }
              }
          }

          // Check for equal max score. If we set a new high score this will set the count to 1
          tmp = cmp_eq(S_curr[col], max_score);
          if (reduce_or(tmp)) {
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set
                  if (extract(i, tmp)) {
                      ++max_count[i];
                      max_pos[i] = closer<uint32_t>(curr, max_pos[i], targets[i]);
                  }
              }
          }

          // new second best score
          tmp = S_curr[col] > sub_score;
          if (reduce_or(tmp)) {
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set, and if far enough
                  if (extract(i, tmp) && max_pos[i] < curr - read_len) {
                      //TODO add -(max_score/sub_score) term
                      insert(extract(i, S_curr[col]), i, sub_score);
                      sub_pos[i] = curr;
                      sub_count[i] = 0;
                  }
              }
          }

          // Repeat sub score
          tmp = cmp_eq(S_curr[col], sub_score);
          if (reduce_or(tmp)) {
              for (uchar i = 0; i < num_reads; ++i) {
                  // Check if the i'th elements MSB is set
                  if (extract(i, tmp)) {
                      ++sub_count[i];
                      sub_pos[i] = closer<uint32_t>(curr, sub_pos[i], targets[i]);
                  }
              }
          }
      }

      /**
       * @brief
       * Extract the i'th element from a vector. No range checking is done.
       * @param i index of element
       * @param vec vector to extract from
       */
      __INLINE__
      NativeT extract(uint8_t i, const CellType<num_reads> &vec) {
          static NativeT *ptr;
          ptr = (NativeT *) &vec;
          return ptr[i];
      }

      /**
       * @brief
       * Insert into the i'th element from a vector. No range checking is done.
       * @param ins element to insert
       * @param i index of element
       * @param vec vector to insert in
       */
      __INLINE__
      void insert(NativeT ins, uint8_t i, const CellType<num_reads> &vec) {
          static NativeT *ptr;
          ptr = (NativeT *) &vec;
          ptr[i] = ins;
      }

    private:

      int read_len; /**< Maximum read length. */

      // Zero vector
      const CellType<num_reads> ZERO_CT = simdpp::splat(0);

      uint8_t
          _match = 2,       /**< Match score, is added */
          _mismatch = 2,    /**< mismatch penalty, is subtracted */
          _gap_open = 3,    /**< gap open penalty, subtracted */
          _gap_extend = 1;  /**< gap extension penalty, subtracted */

      /**
       * Each vector has an 'a' and a 'b' version. Through each row of the
       * matrix fill, their roles are swapped such that one becomes the previous
       * loops data, and the other is filled in.
       * S and D are padded 1 to provide a left column buffer.
       */
      VecType
          *Sa = nullptr, /**< Matrix row */
          *Sb = nullptr,
          *Da = nullptr, /**< Deletion vector */
          *Db = nullptr,
          *Ia = nullptr, /**< Insertion vector */
          *Ib = nullptr;

      CellType<num_reads>
          *S_prev = nullptr,/**< S_prev[n] => S(i-1, n) */
          *S_curr = nullptr,
          *D_prev = nullptr, /**< D_prev[n] => D(i-1, n) */
          *D_curr = nullptr,
          *I_prev = nullptr, /**< I_prev[r] => I(r, j-1) */
          *I_curr = nullptr,
          *swp_tmp;

      CellType<num_reads>
          Ceq, /**< Match score when read_base == ref_base */
          Cneq, /**< mismatch penalty */
          tmp; /**< temporary for use within functions */

      // Optimal alignment info
      CellType<num_reads> max_score;
      std::vector<uint32_t> max_pos;
      std::vector<uint32_t> max_count;

      // Suboptimal alignment info
      CellType<num_reads> sub_score;
      std::vector<uint32_t> sub_pos;
      std::vector<uint32_t> sub_count;

      size_t _max_node_len;

  };

  /**
 * @brief Print alignment in CSV form.
 * @details
 * Print the alignment to os. This ordering matches the way the alignment is parsed
 * from a string. \n
 * read_str,opt_score, opt_alignment_end, opt_cout, sub_score, sub_alignment_end, sub_count \n
 * @param os Output stream
 * @param an Alignment output
 */
  inline std::ostream &operator<<(std::ostream &os, const Alignment &a) {
      os << to_csv(a.read)
          << ',' << a.opt_score
          << ',' << a.opt_align_end
          << ',' << a.opt_count
          << ',' << a.sub_score
          << ',' << a.sub_align_end
          << ',' << a.sub_count
          << ',' << (int) a.corflag;
      return os;
  }

}

TEST_CASE ("Alignment") {

    srand(12345);

        SUBCASE("Node Fill scores") {
        std::vector<std::string> reads;
        reads.push_back("CGT"); // Score 6, Match
        reads.push_back("ATAGCCA"); // Score 10, del
        reads.push_back("ATAACGCCA"); // Score 12, Ins
        reads.push_back("TACCCCA"); // Score 10, mismatch

        Vargas::Graph::Node n;
        n.set_endpos(13);
        n.set_seq("ACGTNATACGCCA");

        Vargas::Aligner<>::AlignmentGroup rb(reads, 16);
        Vargas::Aligner<> a(15, 16);
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
        Vargas::Graph::Node::_newID = 0;
        Vargas::Graph g;

        /**   GGG
        *    /   \
        * AAA     TTTA
        *    \   /
        *     CCC(ref)
        */

        {
            Vargas::Graph::Node n;
            n.set_endpos(3);
            n.set_as_ref();
            std::vector<bool> a = {0, 1, 1};
            n.set_population(a);
            n.set_seq("AAA");
            g.add_node(n);
        }

        {
            Vargas::Graph::Node n;
            n.set_endpos(6);
            n.set_as_ref();
            std::vector<bool> a = {0, 0, 1};
            n.set_population(a);
            n.set_af(0.4);
            n.set_seq("CCC");
            g.add_node(n);
        }

        {
            Vargas::Graph::Node n;
            n.set_endpos(6);
            n.set_not_ref();
            std::vector<bool> a = {0, 1, 0};
            n.set_population(a);
            n.set_af(0.6);
            n.set_seq("GGG");
            g.add_node(n);
        }

        {
            Vargas::Graph::Node n;
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

        std::vector<std::string> reads;
        reads.push_back("CCTT");
        reads.push_back("GGTT");
        reads.push_back("AAGG");
        reads.push_back("AACC");
        reads.push_back("AGGGT");
        reads.push_back("GG");
        reads.push_back("AAATTTA");
        reads.push_back("AAAGCCC");

        Vargas::Aligner<>::AlignmentGroup rb(reads, 8);
        Vargas::Aligner<> a(5, 8);

        Vargas::Aligner<>::Results aligns = a.align(rb, g);
            CHECK(aligns.best_score[0] == 8);
            CHECK(aligns.best_pos[0] == 8);

            CHECK(aligns.best_score[1] == 8);
            CHECK(aligns.best_pos[1] == 8);

            CHECK(aligns.best_score[2] == 8);
            CHECK(aligns.best_pos[2] == 5);

            CHECK(aligns.best_score[3] == 8);
            CHECK(aligns.best_pos[3] == 5);

            CHECK(aligns.best_score[4] == 10);
            CHECK(aligns.best_pos[4] == 7);

            CHECK(aligns.best_score[5] == 4);
            CHECK(aligns.best_pos[5] == 5);

            CHECK(aligns.best_score[6] == 8);
            CHECK(aligns.best_pos[6] == 10);

            CHECK(aligns.best_score[7] == 8);
            CHECK(aligns.best_pos[7] == 6);
    }
}

void node_fill_profile() {
    Vargas::Graph::Node n;
    {
        std::ostringstream ss;
        for (size_t i = 0; i < 10000; ++i) {
            ss << rand_base();
        }
        n.set_seq(ss.str());
    }

    std::vector<std::string> reads;
    for (size_t i = 0; i < 16; ++i) {
        std::ostringstream rd;
        for (size_t r = 0; r < 50; ++r) rd << rand_base();
        reads.push_back(rd.str());
    }

    Vargas::Aligner<>::AlignmentGroup rb(reads, 50);
    Vargas::Aligner<> a(10000, 50);

    std::vector<uint32_t> pos;

    std::cout << "\n10 Node Fill (" << SIMDPP_FAST_INT8_SIZE << "x 50rdlen x 10000bp):\n\t";
    clock_t start = std::clock();
    for (int i = 0; i < 10; ++i) a._test_fill_node(n, rb, pos);
    std::cout << (std::clock() - start) / (double) (CLOCKS_PER_SEC) << " s\n" << std::endl;
}
#endif //VARGAS_ALIGNMENT_H
