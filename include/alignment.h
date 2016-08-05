/**
 * @author Ravi Gaddipati
 * @date June 21, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Aligns groups of short reads to a graph.
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

#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <string>
#include "utils.h"
#include "simdpp/simd.h"
#include "graph.h"
#include "doctest.h"


#define DEBUG_PRINT_SW 0 // Print the full SW matrix for each node aligned
#define DEBUG_PRINT_SW_NUM 0 // Print the matrix for this read number in the alignment group

#define ALIGN_SAM_TYPE_REF "REF"
#define ALIGN_SAM_TYPE_MAXAF "MAXAF"
#define ALIGN_SAM_TYPE_IN "IN"
#define ALIGN_SAM_TYPE_OUT "OUT"
#define ALIGN_SAM_TYPE_TAG "tp"
#define ALIGN_SAM_MAX_POS_TAG "mp"
#define ALIGN_SAM_SUB_POS_TAG "sp"
#define ALIGN_SAM_MAX_SCORE_TAG "ms"
#define ALIGN_SAM_SUB_SCORE_TAG "ss"
#define ALIGN_SAM_MAX_COUNT_TAG "mc"
#define ALIGN_SAM_SUB_COUNT_TAG "sc"
#define ALIGN_SAM_COR_FLAG_TAG "cf"

namespace Vargas {

  /**
   * @brief Main SIMD SW Aligner.
   * @details
   * Aligns a read batch to a reference sequence.
   * Note: "score" means something that is added, "penalty" refers to something
   * that is subtracted. All scores/penalties are provided as positive integers.
   * Most memory is allocated for the class so it can be reused during alignment. To reduce memory usage,
   * the maximum node size can be reduced. \n
   * Usage:\n
   * @code{.cpp}
   * #include "graph.h"
   * #include "alignment.h"
   *
   * Vargas::GraphBuilder gb("reference.fa", "var.bcf");
   * gb.node_len(5);
   * gb.ingroup(100);
   * gb.region("x:0-15");
   *
   * Vargas::Graph g = gb.build();
   * std::vector<std::string> reads = {"ACGT", "GGGG", "ATTA", "CCNT"};
   *
   * Vargas::ByteAligner a(g.max_node_len(), 4);
   * Vargas::ByteAligner::Results res = a.align(reads, g.begin(), g.end());
   *
   * for (int i = 0; i < reads.size(); ++i) {
   *    std::cout << reads[i] << ", score:" << res._max_score[i] << " pos:" << res._max_pos[i] << std::endl;
   * }
   * // Output:
   * // ACGT, score:8 pos:10
   * // GGGG, score:6 pos:65
   * // ATTA, score:8 pos:18
   * // CCNT, score:8 pos:80
   * @endcode
   */

  class ByteAligner {
    public:
      typedef typename std::vector<simdpp::uint8<SIMDPP_FAST_INT8_SIZE>> VecType;

      /**
       * @brief
       * Default constructor uses the following score values: \n
       * Match : 2 \n
       * Mismatch penalty : 2 \n
       * Gap Open penalty : 3 \n
       * Gap Extend penalty : 1 \n
       * @param max_node_len maximum node length
       * @param read_len maximum read length
       */
      ByteAligner(size_t max_node_len,
                  size_t read_len) :
          _read_len(read_len),
          _match_vec(simdpp::splat(2)),
          _mismatch_vec(simdpp::splat(2)),
          _gap_open_extend_vec(simdpp::splat(4)),
          _gap_extend_vec(simdpp::splat(1)),
          _max_node_len(max_node_len),
          _alignment_group(read_len) { _alloc(); }

      /**
       * @brief
       * Set scoring parameters.
       * @param max_node_len max node length
       * @param rlen max read length
       * @param match match score
       * @param mismatch mismatch penalty
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      ByteAligner(size_t max_node_len,
                  size_t read_len,
                  uint8_t match,
                  uint8_t mismatch,
                  uint8_t open,
                  uint8_t extend) :
          _read_len(read_len),
          _match_vec(simdpp::splat(match)),
          _mismatch_vec(simdpp::splat(mismatch)),
          _gap_open_extend_vec(simdpp::splat(open + extend)),
          _gap_extend_vec(simdpp::splat(extend)),
          _max_node_len(max_node_len),
          _alignment_group(read_len) { _alloc(); }

      ~ByteAligner() {
          _dealloc();
      }

      /**
       * @brief
       * Container for a packaged batch of reads.
       * @details
       * Reads are interleaved so each SIMD vector
       * contains bases from all reads, respective to the base number. For example AlignmentGroup[0]
       * would contain the first bases of every read. All reads must be the same length. Minimal error checking.
       * @tparam SIMDPP_FAST_INT8_SIZE max number of reads. If a non-default T is used, this should be set to
       *    SIMDPP_FAST_T_SIZE where T corresponds to the width of T. For ex. Default T=simdpp::uint8 uses
       *    SIMDPP_FAST_INT8_SIZE
       */
      class AlignmentGroup {
        public:

          AlignmentGroup(size_t read_len) : _read_len(read_len), _packaged_reads(read_len) {}

          __INLINE__ void load_reads(const std::vector<std::string> &reads, size_t begin, size_t end) {
              load_reads(std::vector<std::string>(reads.begin() + begin, reads.begin() + end));
          }

          /**
           * @param batch load the given vector of reads.
           */
          __INLINE__ void load_reads(const std::vector<std::string> &batch) {
              std::vector<std::vector<Base>> _reads;
              for (auto &b : batch) _reads.push_back(seq_to_num(b));
              load_reads(_reads);
          }

          /**
           * @param batch load the given vector of reads.
           */
          __INLINE__ void load_reads(const std::vector<std::vector<Base>> &batch) {
              _package_reads(batch);
          }

          /**
           * @brief
           * Return the i'th base of every read in a simdpp vector.
           * @param i base index.
           */
          const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &at(int i) const {
              return _packaged_reads.at(i);
          }

          /**
           * @brief
           * Pointer to raw packaged read data.
           * @return simdpp::uint8<SIMDPP_FAST_INT8_SIZE> pointer
           */
          const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> *data() const {
              return _packaged_reads.data();
          }

          /**
           * @brief
           * Non const version of at(i).
           * @param i base index
           */
          simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &operator[](int i) {
              return _packaged_reads.at(i);
          }

          /**
           * @brief
           * Returns optimal number of reads in a batch based on SIMD architecture.
           * @return batch size.
           */
          size_t group_size() const { return SIMDPP_FAST_INT8_SIZE; }

          /**
           * @return iterator to the beginning of the packaged reads.
           */
          typename std::vector<simdpp::uint8<SIMDPP_FAST_INT8_SIZE>>::const_iterator begin() const {
              return _packaged_reads.begin();
          }

          /**
           * @return iterator to the end of the packaged reads.
           */
          typename std::vector<simdpp::uint8<SIMDPP_FAST_INT8_SIZE>>::const_iterator end() const {
              return _packaged_reads.end();
          }

        private:

          const size_t _read_len;

          /**
           * _packaged_reads[i] contains all i'th bases.
           * The length of _packaged_reads is the length of the read,
           * where as the length of _packaged_reads[i] is the number
           * of reads.
           */
          std::vector<simdpp::uint8<SIMDPP_FAST_INT8_SIZE>> _packaged_reads;

          /**
           * Interleaves reads so all same-index base positions are in one
           * vector. Empty spaces are padded with Base::N.
           * @param _reads vector of reads to package
           */
          __INLINE__ void _package_reads(const std::vector<std::vector<Base>> &_reads) {
              assert(_reads.size() <= SIMDPP_FAST_INT8_SIZE);
              // Interleave reads
              // For each read (read[i] is in _packaged_reads[0..n][i]
              for (size_t r = 0; r < _reads.size(); ++r) {
                  assert(_reads[r].size() == _read_len);
                  // Put each base in the appropriate vector element
                  for (size_t p = 0; p < _read_len; ++p) {
                      insert(_reads[r][p], r, _packaged_reads[p]);
                  }
              }

              // Pad underful batches
              for (size_t r = _reads.size(); r < SIMDPP_FAST_INT8_SIZE; ++r) {
                  for (size_t p = 0; p < _read_len; ++p) {
                      insert(Base::N, r, _packaged_reads[p]);
                  }
              }

          }

      };

      /**
       * @brief
       * Struct to return the alignment results
       */
      struct Results {
          std::vector<uint32_t> max_pos;
          /**< Best positions */
          std::vector<uint32_t> sub_pos;
          /**< Second best positions */

          std::vector<uint8_t> max_count;
          /**< Occurances of max_pos */
          std::vector<uint8_t> sub_count;
          /**< Occurances of _sub_pos */

          std::vector<uint8_t> max_score;
          /**< Best scores */
          std::vector<uint8_t> sub_score;
          /**< Second best scores */

          std::vector<uint8_t> cor_flag;
          /**< 1 for target matching best score, 2 for matching sub score, 0 otherwise */

          /**
           * @brief
           * Resize all result vectors.
           */
          void resize(size_t size) {
              max_pos.resize(size);
              sub_pos.resize(size);
              max_count.resize(size);
              sub_count.resize(size);
              max_score.resize(size);
              sub_score.resize(size);
              cor_flag.resize(size);
          }
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
          using namespace simdpp;

          _match_vec = splat(match);
          _mismatch_vec = splat(mismatch);
          _gap_open_extend_vec = splat(open + extend);
          _gap_extend_vec = splat(extend);

      }

      /**
       * @return maximum number of reads that can be aligned at once.
       */
      inline size_t read_capacity() const { return SIMDPP_FAST_INT8_SIZE; }

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group vector of reads to align to
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @return Results packet
       */
      Results align(const std::vector<std::string> &read_group,
                    Graph::FilteringIter begin,
                    Graph::FilteringIter end) {
          std::vector<uint32_t> targets(read_group.size());
          std::fill(targets.begin(), targets.end(), 0);
          return align(read_group, targets, begin, end);
      }

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group vector of reads to align to
       * @param targets Keep cell scores for the given positions
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @return Results packet
       */
      Results align(const std::vector<std::string> &read_group,
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
       * @param read_group vector of reads to align to
       * @param targets Origin of read, determines cor_flag
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @param aligns Results packet to populate
       */
      inline void align_into(const std::vector<std::string> read_group,
                             std::vector<uint32_t> targets,
                             Graph::FilteringIter begin,
                             Graph::FilteringIter end,
                             Results &aligns) {
          using namespace simdpp;

          if (targets.size() != read_group.size()) {
              std::cerr << "Warning: target size does not match read group size, setting to zero." << std::endl;
              targets.resize(read_group.size());
              std::fill(targets.begin(), targets.end(), 0);
          }

          size_t num_groups = read_group.size() / SIMDPP_FAST_INT8_SIZE; // Full alignment groups
          aligns.resize(read_group.size());
          std::fill(aligns.cor_flag.begin(), aligns.cor_flag.end(), 0);

          std::unordered_map<uint32_t, _seed> seed_map; // Maps node ID's to the ending matrix columns of the node
          _seed seed(_read_len), nxt(_read_len);

          for (size_t group = 0; group < num_groups; ++group) {
              seed_map.clear();

              // Subset of read set
              const size_t beg_offset = group * SIMDPP_FAST_INT8_SIZE;
              const size_t end_offset = (group + 1) * SIMDPP_FAST_INT8_SIZE;
              _alignment_group.load_reads(read_group, beg_offset, end_offset);
              _max_score = ZERO_CT;
              _sub_score = ZERO_CT;

              _max_pos = aligns.max_pos.data() + beg_offset;
              _sub_pos = aligns.sub_pos.data() + beg_offset;
              _max_count = aligns.max_count.data() + beg_offset;
              _sub_count = aligns.sub_count.data() + beg_offset;
              _cor_flag = aligns.cor_flag.data() + beg_offset;
              _targets = targets.data() + beg_offset;

              for (auto gi = begin; gi != end; ++gi) {
                  _get_seed(gi.incoming(), seed_map, &seed);
                  if (gi->is_pinched()) seed_map.clear();
                  _fill_node(*gi, _alignment_group, &seed, &nxt);
                  seed_map.emplace(gi->id(), nxt);
              }

              memcpy(aligns.max_score.data() + beg_offset, &_max_score, SIMDPP_FAST_INT8_SIZE * sizeof(uint8_t));
              memcpy(aligns.sub_score.data() + beg_offset, &_sub_score, SIMDPP_FAST_INT8_SIZE * sizeof(uint8_t));
          }

          // If we need a padded group at the end
          if (read_group.size() % SIMDPP_FAST_INT8_SIZE) {
              seed_map.clear();

              size_t len = read_group.size() - (num_groups * SIMDPP_FAST_INT8_SIZE);
              size_t offset = num_groups * SIMDPP_FAST_INT8_SIZE;

              _alignment_group.load_reads(read_group, num_groups * SIMDPP_FAST_INT8_SIZE, read_group.size());
              _max_score = ZERO_CT;
              _sub_score = ZERO_CT;

              std::vector<uint32_t>
                  tmp_targets(SIMDPP_FAST_INT8_SIZE),
                  tmp_max_pos(SIMDPP_FAST_INT8_SIZE),
                  tmp_sub_pos(SIMDPP_FAST_INT8_SIZE);
              std::vector<uint8_t>
                  tmp_max_count(SIMDPP_FAST_INT8_SIZE),
                  tmp_sub_count(SIMDPP_FAST_INT8_SIZE),
                  tmp_cor_flag(SIMDPP_FAST_INT8_SIZE);

              _max_pos = tmp_max_pos.data();
              _sub_pos = tmp_sub_pos.data();
              _max_count = tmp_max_count.data();
              _sub_count = tmp_sub_count.data();
              _cor_flag = tmp_cor_flag.data();
              _targets = tmp_targets.data();

              memcpy(_targets, targets.data() + offset, len * sizeof(uint32_t));


              for (auto &gi = begin; gi != end; ++gi) {
                  _get_seed(gi.incoming(), seed_map, &seed);
                  if (gi->is_pinched()) seed_map.clear();
                  _fill_node(*gi, _alignment_group, &seed, &nxt);
                  seed_map.emplace(gi->id(), nxt);
              }

              memcpy(aligns.max_score.data() + offset, &_max_score, len * sizeof(uint8_t));
              memcpy(aligns.sub_score.data() + offset, &_sub_score, len * sizeof(uint8_t));

              memcpy(aligns.max_pos.data() + offset, tmp_max_pos.data(), len * sizeof(uint32_t));
              memcpy(aligns.max_count.data() + offset, tmp_max_count.data(), len * sizeof(uint8_t));
              memcpy(aligns.sub_pos.data() + offset, tmp_sub_pos.data(), len * sizeof(uint32_t));
              memcpy(aligns.sub_count.data() + offset, tmp_sub_count.data(), len * sizeof(uint8_t));
              memcpy(aligns.cor_flag.data() + offset, tmp_cor_flag.data(), len * sizeof(uint8_t));
          }

      }

      /**
       * @brief
       * Ensures that the graph is topographically sorted.
       * @param begin Graph begin iterator
       * @param end Graph end iterator
       */
      bool validate(Graph::FilteringIter begin,
                    Graph::FilteringIter end) {
          std::unordered_set<size_t> filled;
          bool ret = true;
          for (auto &gi = begin; gi != end; ++gi) {
              filled.insert(gi->id());
              for (auto i : gi.incoming()) {
                  if (filled.count(i) == 0) {
                      ret = false;
                      std::cerr << "Node (ID:" << gi->id() << ", POS:" << gi->end() << ")"
                                << " hit before previous node " << i << std::endl;
                  }
              }
          }
          return ret;
      }

    private:
      /**
       * @brief
       * Ending vectors from a previous node
       */
      struct _seed {
          _seed(int _read_len) : S_col(_read_len), I_col(_read_len) { }
          VecType S_col;
          /**< Last column of score matrix.*/
          VecType I_col;
          /**< Last column of I vector.*/
      };

      /**
       * @brief
       * Returns the best seed from all previous nodes.
       * @param prev_ids All nodes preceding _curr_posent node
       * @param seed_map ID->seed map for all previous nodes
       * @param seed best seed to populate
       * @throws std::logic_error if a node listed as a previous node but it has not been encountered yet. i.e. not topographically sorted.
       */
      __INLINE__
      void _get_seed(const std::vector<uint32_t> &prev_ids,
                     const std::unordered_map<uint32_t, _seed> &seed_map,
                     _seed *seed) {
          using namespace simdpp;

          const _seed *ns;

          try {
              for (size_t i = 0; i < _read_len; ++i) {
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
          if (_Sa) delete _Sa;
          if (_Sb) delete _Sb;
          if (_Da) delete _Da;
          if (_Db) delete _Db;
          if (_Ia) delete _Ia;
          if (_Ib) delete _Ib;

          _Sa = nullptr;
          _Sb = nullptr;
          _Da = nullptr;
          _Db = nullptr;
          _Ia = nullptr;
          _Ib = nullptr;

          _S_prev = nullptr;
          _S_curr = nullptr;
          _D_prev = nullptr;
          _D_curr = nullptr;
          _I_prev = nullptr;
          _I_curr = nullptr;
      }

      /**
       * @brief
       * Allocate S and D vectors.
       */
      void _alloc() {
          _Sa = new VecType(_max_node_len);
          _Sb = new VecType(_max_node_len);
          _Da = new VecType(_max_node_len);
          _Db = new VecType(_max_node_len);
          _Ia = new VecType(_read_len);
          _Ib = new VecType(_read_len);

          _S_prev = _Sa->data();
          _S_curr = _Sb->data();
          _D_prev = _Da->data();
          _D_curr = _Db->data();
          _I_prev = _Ia->data();
          _I_curr = _Ib->data();
      }

      /**
       * @brief
       * Computes local alignment to the node.
       * @param n Node to align to
       * @param read_group AlignmentGroup to align
       * @param s seeds from previous nodes
       * @param nxt seed for next nodes
       */
      void _fill_node(const Graph::Node &n,
                      const AlignmentGroup &read_group,
                      const _seed *s,
                      _seed *nxt) {

          assert(n.seq().size() <= _max_node_len);

          const Base *node_seq = n.seq().data();
          const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> *read_ptr = read_group.data();
          const size_t seq_size = n.seq().size();
          const size_t node_origin = n.end() - seq_size + 2;

          // top left corner
          _fill_cell_rzcz(read_ptr[0], node_seq[0], s);
          _fill_cell_finish(0, 0, node_origin);

          // top row
          for (uint32_t c = 1; c < seq_size; ++c) {
              _fill_cell_rz(read_ptr[0], node_seq[c], c);
              _fill_cell_finish(0, c, node_origin);
          }

          nxt->S_col[0] = _S_curr[seq_size - 1];

          // Rest of the rows
          for (uint32_t r = 1; r < _read_len; ++r) {
              // Swap the rows we are filling in. The previous row/col becomes what we fill in.
              _swp_tmp0 = _S_prev;
              _S_prev = _S_curr;
              _S_curr = _swp_tmp0;

              _swp_tmp0 = _D_prev;
              _D_prev = _D_curr;
              _D_curr = _swp_tmp0;

              _swp_tmp0 = _I_prev;
              _I_prev = _I_curr;
              _I_curr = _swp_tmp0;

              // first col
              _fill_cell_cz(read_ptr[r], node_seq[0], r, s);
              _fill_cell_finish(r, 0, node_origin);

              // Inner grid
              for (uint32_t c = 1; c < seq_size; ++c) {
                  _fill_cell(read_ptr[r], node_seq[c], r, c);
                  _fill_cell_finish(r, c, node_origin);
              }

              nxt->S_col[r] = _S_curr[seq_size - 1];

          }

          // origin vector of what is now _I_curr
          nxt->I_col = _Ia->data() == _I_curr ? *_Ia : *_Ib;
      }

      /**
       * @brief
       * Fills the top left cell.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param s alignment seed from previous node
       */
      __INLINE__
      void _fill_cell_rzcz(const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &read_base,
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
       * @param col _curr_posent column in matrix
       */
      __INLINE__
      void _fill_cell_rz(const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &read_base,
                         const Base &ref,
                         const uint32_t &col) {
          _D(col, ZERO_CT, ZERO_CT);
          _I(0, _S_curr[col - 1]);
          _M(col, read_base, ref, ZERO_CT);
      }

      /**
       * @brief
       * Fills cells when col is 0.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row _curr_posent row in matrix
       * @param s alignment seed from previous node
       */
      __INLINE__
      void _fill_cell_cz(const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &read_base,
                         const Base &ref,
                         const uint32_t &row,
                         const _seed *s) {
          _D(0, _D_prev[0], _S_prev[0]);
          _I(row, s->S_col[row]);
          _M(0, read_base, ref, s->S_col[row - 1]);
      }

      /**
       * @brief
       * Fills the _curr_posent cell.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row _curr_posent row in matrix
       * @param col _curr_posent column in matrix
       */
      __INLINE__
      void _fill_cell(const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &read_base,
                      const Base &ref,
                      const uint32_t &row,
                      const uint32_t &col) {

          _D(col, _D_prev[col], _S_prev[col]);
          _I(row, _S_curr[col - 1]);
          _M(col, read_base, ref, _S_prev[col - 1]);
      }

      /**
       * @brief
       * Score if there is a deletion
       * @param col _curr_posent column
       * @param Dp Previous D value at _curr_posent col.
       * @param Sp Previous S value at _curr_posent col.
       */
      __INLINE__
      void _D(const uint32_t &col,
              const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &Dp,
              const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &Sp) {
          using namespace simdpp;

          // D(i,j) = D(i-1,j) - gap_extend
          // Dp is _D_prev[col], 0 for row=0
          _D_curr[col] = sub_sat(Dp, _gap_extend_vec);   // _tmp0 = S(i-1,j) - ( gap_open + gap_extend)
          // Sp is _S_prev[col], 0 for row=0
          // D(i,j) = max{ D(i-1,j) - gap_extend, S(i-1,j) - ( gap_open + gap_extend) }
          _tmp0 = sub_sat(Sp, _gap_open_extend_vec);
          _D_curr[col] = max(_D_curr[col], _tmp0);

      }

      /**
       * @brief
       * Score if there is an insertion
       * @param row _curr_posent row
       * @param Sc Previous S value (cell to the left)
       */
      __INLINE__
      void _I(const uint32_t &row,
              const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &Sc) {
          using namespace simdpp;

          // I(i,j) = I(i,j-1) - gap_extend
          _I_curr[row] = sub_sat(_I_prev[row], _gap_extend_vec);  // I: I(i,j-1) - gap_extend
          // _tmp0 = S(i,j-1) - (gap_open + gap_extend)
          // Sc is _S_curr[col - 1], seed->S_col[row] for col=0
          _tmp0 = sub_sat(Sc, _gap_open_extend_vec);
          _I_curr[row] = max(_I_curr[row], _tmp0);

      }

      /**
       * @brief
       * Best score if there is a match/mismatch. Uses _S_prev.
       * @param col _curr_posent column
       * @param read read base vector
       * @param ref reference sequence base
       * @param Sp Previous S val at col-1 (upper left cell)
       */
      __INLINE__
      void _M(uint32_t col,
              const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &read,
              const Base &ref,
              const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &Sp) {
          using namespace simdpp;

          if (ref != Base::N) {
              _Ceq = ZERO_CT;
              _Cneq = ZERO_CT;

              // Set all mismatching pairs to _mismatch
              _tmp0 = cmp_neq(read, ref);
              _Cneq = _tmp0 & _mismatch_vec;   // If the read base is Base::N, set to 0 (_Ceq)
              _tmp0 = cmp_eq(read, _N_VEC);
              _Cneq = blend(_Ceq, _Cneq, _tmp0);

              // b is not N, so all equal bases are valid
              _tmp0 = cmp_eq(read, ref);
              _Ceq = _tmp0 & _match_vec;

              // Sp is _S_prev[col - 1], 0 for row=0
              // Seed->S_col[row - 1] for col=0
              _S_curr[col] = add_sat(Sp, _Ceq);   // Add match scores
              _S_curr[col] = sub_sat(_S_curr[col], _Cneq); // Subtract mismatch scores
          }

      }


      /**
       * @brief
       * Takes the max of D,I, and M vectors and stores the _curr_posent best score/position
       * Currently does not support non-default template args
       * @param row _curr_posent row
       * @param col _curr_posent column
       * @param node_origin Current position, used to get absolute alignment position
       */
      __INLINE__ __UNROLL__
      void _fill_cell_finish(const uint32_t &row,
                             const uint32_t &col,
                             const uint32_t &node_origin) {
          using namespace simdpp;



          // S(i,j) = max{ D(i,j), I(i,j), S(i-1,j-1) + C(s,t) }
          _S_curr[col] = max(_D_curr[col], _S_curr[col]);
          _S_curr[col] = max(_I_curr[row], _S_curr[col]);

          _curr_pos = node_origin + col;    // absolute position in reference sequence

          _tmp0 = _S_curr[col] > _max_score;
          if (extract_bits_any(_tmp0)) {
              // Check for new or equal high scores
              _max_score = max(_S_curr[col], _max_score);
              for (int i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  if (_tmp0_ptr[i]) {
                      // Demote old max to submax
                      if (_curr_pos > _max_pos[i] + _read_len) {
                          _sub_score_ptr[i] = _max_score_ptr[i];
                          _sub_pos[i] = _max_pos[i];
                          _sub_count[i] = _max_count[i];
                          _cor_flag[i] = (_cor_flag[i] == 1) * 2;
                      }
                      _max_count[i] = 1;
                      _cor_flag[i] = (_cor_flag[i] == 2) * 2;
                  }
              }
          }

          _tmp0 = cmp_eq(_S_curr[col], _max_score);
          if (extract_bits_any(_tmp0)) {
              // Check for equal max score.
              for (uint8_t i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  if (_tmp0_ptr[i]) {
                      _max_count[i] += _curr_pos > (_max_pos[i] + _read_len);
                      _max_pos[i] = _curr_pos;
                      if (_max_pos[i] == _targets[i]) _cor_flag[i] = 1;
                  }
              }
          }


          // Greater than old sub max and less than max score (prevent repeats of max triggering)
          _tmp0 = (_S_curr[col] > _sub_score) & (_S_curr[col] < _max_score);
          if (extract_bits_any(_tmp0)) {
              // new second best score
              for (uint8_t i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  if (_tmp0_ptr[i] && _curr_pos > _max_pos[i] + _read_len) {
                      _sub_score_ptr[i] = extract(i, _S_curr[col]);
                      _sub_count[i] = 1;
                      _sub_pos[i] = _curr_pos;
                      if (_curr_pos == _targets[i]) _cor_flag[i] = 2;
                      else _cor_flag[i] = _cor_flag[i] == 1;
                  }
              }
          }

          _tmp0 = cmp_eq(_S_curr[col], _sub_score);
          if (extract_bits_any(_tmp0)) {
              // Repeat sub score
              for (uint8_t i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  if (_tmp0_ptr[i] && _curr_pos > _max_pos[i] + _read_len) {
                      _sub_count[i] += _curr_pos > (_sub_pos[i] + _read_len);
                      _sub_pos[i] = _curr_pos;
                      if (_curr_pos == _targets[i]) _cor_flag[i] = 2;
                  }
              }
          }
      }

      /*********************************** Variables ***********************************/

      size_t _read_len; /**< Maximum read length. */

      // Zero vector
      const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> ZERO_CT = simdpp::splat(0);
      const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> _N_VEC = simdpp::splat(Base::N);

      simdpp::uint8<SIMDPP_FAST_INT8_SIZE>
          _match_vec,
          _mismatch_vec,
          _gap_open_extend_vec,
          _gap_extend_vec;

      /**
       * Each vector has an 'a' and a 'b' version. Through each row of the
       * matrix fill, their roles are swapped such that one becomes the previous
       * loops data, and the other is filled in.
       * S and D are padded 1 to provide a left column buffer.
       */
      VecType
          *_Sa = nullptr,    /**< Matrix row */
          *_Sb = nullptr,
          *_Da = nullptr,    /**< Deletion vector */
          *_Db = nullptr,
          *_Ia = nullptr,    /**< Insertion vector */
          *_Ib = nullptr;

      simdpp::uint8<SIMDPP_FAST_INT8_SIZE>
          *_S_prev = nullptr,    /**< _S_prev[n] => S(i-1, n) */
          *_S_curr = nullptr,
          *_D_prev = nullptr,    /**< _D_prev[n] => D(i-1, n) */
          *_D_curr = nullptr,
          *_I_prev = nullptr,    /**< _I_prev[r] => I(r, j-1) */
          *_I_curr = nullptr,
          *_swp_tmp0;

      simdpp::uint8<SIMDPP_FAST_INT8_SIZE>
          _tmp0, /**< temporary for use within functions */
          _Ceq,  /**< Match score when read_base == ref_base */
          _Cneq; /**< mismatch penalty */

      uint32_t _curr_pos;

      // Optimal alignment info
      simdpp::uint8<SIMDPP_FAST_INT8_SIZE> _max_score;
      uint32_t *_max_pos;
      uint8_t *_max_count;

      // Suboptimal alignment info
      simdpp::uint8<SIMDPP_FAST_INT8_SIZE> _sub_score;
      uint32_t *_sub_pos;
      uint8_t *_sub_count;

      uint8_t *_tmp0_ptr = (uint8_t *) &_tmp0;
      uint8_t *_max_score_ptr = (uint8_t *) &_max_score;
      uint8_t *_sub_score_ptr = (uint8_t *) &_sub_score;

      uint8_t *_cor_flag;
      uint32_t *_targets;

      size_t _max_node_len;

      AlignmentGroup _alignment_group;

  };

}

TEST_CASE ("Alignment") {

        SUBCASE("Graph Alignment") {
        Vargas::Graph::Node::_newID = 0;
        Vargas::Graph g;

        /**
        *     GGG
        *    /   \
        * AAA     TTTA
        *    \   /
        *     CCC(ref)
        */

        {
            Vargas::Graph::Node n;
            n.set_endpos(2);
            n.set_as_ref();
            std::vector<bool> a = {0, 1, 1};
            n.set_population(a);
            n.set_seq("AAA");
            g.add_node(n);
        }

        {
            Vargas::Graph::Node n;
            n.set_endpos(5);
            n.set_as_ref();
            std::vector<bool> a = {0, 0, 1};
            n.set_population(a);
            n.set_af(0.4);
            n.set_seq("CCC");
            g.add_node(n);
        }

        {
            Vargas::Graph::Node n;
            n.set_endpos(5);
            n.set_not_ref();
            std::vector<bool> a = {0, 1, 0};
            n.set_population(a);
            n.set_af(0.6);
            n.set_seq("GGG");
            g.add_node(n);
        }

        {
            Vargas::Graph::Node n;
            n.set_endpos(9);
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
        reads.push_back("NNNCCTT");
        reads.push_back("NNNGGTT");
        reads.push_back("NNNAAGG");
        reads.push_back("NNNAACC");
        reads.push_back("NNAGGGT");
        reads.push_back("NNNNNGG");
        reads.push_back("AAATTTA");
        reads.push_back("AAAGCCC");

        Vargas::ByteAligner a(5, 7);

        std::vector<uint32_t> origins = {8, 8, 5, 5, 7, 6, 10, 4};
        Vargas::ByteAligner::Results aligns = a.align(reads, origins, g.begin(), g.end());
            CHECK(aligns.max_score[0] == 8);
            CHECK(aligns.max_pos[0] == 8);
            CHECK((int) aligns.cor_flag[0] == 1);

            CHECK(aligns.max_score[1] == 8);
            CHECK(aligns.max_pos[1] == 8);
            CHECK((int) aligns.cor_flag[1] == 1);

            CHECK(aligns.max_score[2] == 8);
            CHECK(aligns.max_pos[2] == 5);
            CHECK((int) aligns.cor_flag[2] == 1);

            CHECK(aligns.max_score[3] == 8);
            CHECK(aligns.max_pos[3] == 5);
            CHECK((int) aligns.cor_flag[3] == 1);

            CHECK(aligns.max_score[4] == 10);
            CHECK(aligns.max_pos[4] == 7);
            CHECK((int) aligns.cor_flag[4] == 1);

            CHECK(aligns.max_score[5] == 4);
            CHECK(aligns.max_pos[5] == 6);
            CHECK((int) aligns.cor_flag[5] == 1);

            CHECK(aligns.max_score[6] == 8);
            CHECK(aligns.max_pos[6] == 10);
            CHECK((int) aligns.cor_flag[6] == 1);

            CHECK(aligns.max_score[7] == 8);
            CHECK(aligns.max_pos[7] == 4);
            CHECK((int) aligns.cor_flag[7] == 1);

    }
}
#endif //VARGAS_ALIGNMENT_H
