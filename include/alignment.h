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

#define ALIGN_ENABLE_PROFILE 0

#if ALIGN_ENABLE_PROFILE

#include <unordered_map>
#include <ctime>
std::unordered_map<std::string, time_t> _PROF_TOTALS;
std::unordered_map<std::string, time_t> _PROF_STARTS;

#define PROF_BEGIN(name) { _PROF_STARTS[name] = std::clock();}
#define PROF_END(name) {_PROF_TOTALS[name] += std::clock() - _PROF_STARTS[name];}
#define PROF_TERMINATE { \
                        time_t TOT = 0; \
                        for (auto &p : _PROF_TOTALS) TOT += p.second; \
                        for (auto &p : _PROF_TOTALS) { \
                            std::cerr << p.first << ": " << ((double) p.second / (CLOCKS_PER_SEC/1000)) << "ms, " \
                            << ((double) (100 * p.second) / TOT) << "%" << std::endl; \
                        } \
                        _PROF_TOTALS.clear(); \
                        _PROF_STARTS.clear(); \
                        }

#else

#define PROF_BEGIN(name)
#define PROF_END(name)
#define PROF_TERMINATE

#endif

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
       * Mismatch : -2 \n
       * Gap Open : 3 \n
       * Gap Extend : 1 \n
       * @param max_node_len maximum node length
       * @param read_len maximum read length
       */
      ByteAligner(size_t max_node_len, size_t read_len) :
          _read_len(read_len),
          _max_node_len(max_node_len) { _alloc(); }

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
          _match(match), _mismatch(mismatch), _gap_open(open), _gap_extend(extend),
          _max_node_len(max_node_len) { _alloc(); }

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

          AlignmentGroup() { }

          /**
           * @brief
           * Read length is set to first read size.
           * @param batch package the given vector of reads. Must be nonempty.
           */
          AlignmentGroup(const std::vector<std::vector<Base>> &batch) {
              load_reads(batch);
          }

          /**
           * @brief
           * Read length is set to first read size.
           * @param batch package the given vector of reads. Must be nonempty.
           */
          AlignmentGroup(const std::vector<std::string> &batch) {
              load_reads(batch);
          }

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
              _read_len = batch[0].size();
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

          size_t _read_len;

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
              PROF_BEGIN("Package Reads")
              _packaged_reads.resize(_read_len);

              // allocate memory
              uint8_t **pckg = (uint8_t **) malloc(_read_len * sizeof(uint8_t *));
              for (size_t i = 0; i < _read_len; ++i) {
                  pckg[i] = (uint8_t *) malloc(SIMDPP_FAST_INT8_SIZE * sizeof(uint8_t));
              }

              // Interleave reads
              // For each read (read[i] is in _packaged_reads[0..n][i]
              for (size_t r = 0; r < _reads.size(); ++r) {
                  assert(_reads[r].size() == _read_len);
                  // Put each base in the appropriate vector element
                  for (size_t p = 0; p < _read_len; ++p) {
                      pckg[p][r] = _reads[r][p];
                  }
              }

              // Pad underful batches
              for (size_t r = _reads.size(); r < SIMDPP_FAST_INT8_SIZE; ++r) {
                  for (size_t p = 0; p < _read_len; ++p) {
                      pckg[p][r] = Base::N;
                  }
              }

              // Load into vectors
              for (size_t i = 0; i < _read_len; ++i) {
                  _packaged_reads[i] = simdpp::load(pckg[i]);
                  free(pckg[i]);
              }
              free(pckg);

              PROF_END("Package Reads")
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
          _match = match;
          _mismatch = mismatch;
          _gap_open = open;
          _gap_extend = extend;
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

          if (!validate(begin, end)) {
              throw std::logic_error("Graph not topographically ordered. Aborting.");
          }

          size_t num_groups = read_group.size() / SIMDPP_FAST_INT8_SIZE;
          aligns.resize(read_group.size());

          std::unordered_map<uint32_t, _seed> seed_map; // Maps node ID's to the ending columns of the node

          for (size_t group = 0; group < num_groups; ++group) {
              seed_map.clear();
              _seed seed(_read_len);

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

              for (auto gi = begin; gi != end; ++gi) {
                  _get_seed(gi.incoming(), seed_map, &seed);
                  if (gi->is_pinched()) seed_map.clear();
                  seed_map.emplace(gi->id(), _fill_node(*gi, _alignment_group, &seed, targets));
              }

              memcpy(aligns.max_score.data() + beg_offset, &_max_score, SIMDPP_FAST_INT8_SIZE * sizeof(uint8_t));
              memcpy(aligns.sub_score.data() + beg_offset, &_sub_score, SIMDPP_FAST_INT8_SIZE * sizeof(uint8_t));
          }

          // If we need a padded group at the end
          if (read_group.size() % SIMDPP_FAST_INT8_SIZE) {
              seed_map.clear();
              _seed seed(_read_len);

              _alignment_group.load_reads(read_group, num_groups * SIMDPP_FAST_INT8_SIZE, read_group.size());
              _max_score = ZERO_CT;
              _sub_score = ZERO_CT;

              std::vector<uint32_t>
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

              for (auto &gi = begin; gi != end; ++gi) {
                  _get_seed(gi.incoming(), seed_map, &seed);
                  if (gi->is_pinched()) seed_map.clear();
                  seed_map.emplace(gi->id(), _fill_node(*gi, _alignment_group, &seed, targets));
              }

              size_t len = read_group.size() - (num_groups * SIMDPP_FAST_INT8_SIZE);
              size_t offset = num_groups * SIMDPP_FAST_INT8_SIZE;
              memcpy(aligns.max_score.data() + offset, &_max_score, len * sizeof(uint8_t));
              memcpy(aligns.sub_score.data() + offset, &_sub_score, len * sizeof(uint8_t));

              memcpy(aligns.max_pos.data() + offset, tmp_max_pos.data(), len * sizeof(uint32_t));
              memcpy(aligns.max_count.data() + offset, tmp_max_count.data(), len * sizeof(uint8_t));
              memcpy(aligns.sub_pos.data() + offset, tmp_sub_pos.data(), len * sizeof(uint32_t));
              memcpy(aligns.sub_count.data() + offset, tmp_sub_count.data(), len * sizeof(uint8_t));
              memcpy(aligns.cor_flag.data() + offset, tmp_cor_flag.data(), len * sizeof(uint8_t));
          }

          PROF_TERMINATE
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
       */
      void _get_seed(const std::vector<uint32_t> &prev_ids,
                     const std::unordered_map<uint32_t, _seed> &seed_map,
                     _seed *seed) {
          PROF_BEGIN("Get Seed")
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
          PROF_END("Get Seed")
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
       * Computes local alignment of the node, with no previous seed.
       * @param n Node to align to
       * @param read_group AlignmentGroup to align
       */
      _seed _fill_node(const Graph::Node &n,
                       const AlignmentGroup &read_group,
                       const std::vector<uint32_t> &targets) {
          _seed s(_read_len);
          s = _fill_node(n, read_group, &s, targets);
          return s;
      }

      /**
       * @brief
       * Computes local alignment to the node.
       * @param n Node to align to
       * @param read_group AlignmentGroup to align
       * @param s seeds from previous nodes
       */
      _seed _fill_node(const Graph::Node &n,
                       const AlignmentGroup &read_group,
                       const _seed *s,
                       const std::vector<uint32_t> &targets) {

          assert(n.seq().size() < _max_node_len);

          _seed nxt(_read_len);  // Seed for next node
          const Base *node_seq = n.seq().data();
          const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> *read_ptr = read_group.data();
          const size_t seq_size = n.seq().size();
          const size_t node_origin = n.end() - seq_size;

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
                  std::cout << (int) simdpp::extract<DEBUG_PRINT_SW_NUM>(_S_curr[i]) << '\t';
              std::cout << "MAX| ";
              for (auto p : _max_pos) std::cout << p << '\t';
              std::cout << "SUB| ";
              for (auto p : _sub_pos) std::cout << p << '\t';
              std::cout << std::endl;
          #endif

          nxt.S_col[0] = _S_curr[seq_size - 1];

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
              _fill_cell_finish(r, 0, node_origin, targets);

              // Inner grid
              for (uint32_t c = 1; c < seq_size; ++c) {
                  _fill_cell(read_ptr[r], node_seq[c], r, c);
                  _fill_cell_finish(r, c, node_origin, targets);
              }

              nxt.S_col[r] = _S_curr[seq_size - 1];

              #if DEBUG_PRINT_SW
              std::cout << num_to_base((Base) simdpp::extract<DEBUG_PRINT_SW_NUM>(read_ptr[r])) << '\t';
                  for (size_t i = 0; i < n.seq().size(); ++i)
                      std::cout << (int) simdpp::extract<DEBUG_PRINT_SW_NUM>(_S_curr[i]) << '\t';
                  std::cout << "MAX| ";
                  for (auto p : _max_pos) std::cout << p << '\t';
                  std::cout << "SUB| ";
                  for (auto p : _sub_pos) std::cout << p << '\t';
                  std::cout << std::endl;
              #endif

          }

          // origin vector of what is now _I_curr
          nxt.I_col = (_Ia->data() == _I_curr) ? *_Ia : *_Ib;
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
          using namespace simdpp;

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
          PROF_BEGIN("D")

          // D(i,j) = D(i-1,j) - gap_extend
          // Dp is _D_prev[col], 0 for row=0
          _D_curr[col] = sub_sat(Dp, _gap_extend);   // _tmp0 = S(i-1,j) - ( gap_open + gap_extend)
          // Sp is _S_prev[col], 0 for row=0
          // D(i,j) = max{ D(i-1,j) - gap_extend, S(i-1,j) - ( gap_open + gap_extend) }
          _tmp0 = sub_sat(Sp, _gap_extend + _gap_open);
          _D_curr[col] = max(_D_curr[col], _tmp0);

          PROF_END("D")
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
          PROF_BEGIN("I")

          // I(i,j) = I(i,j-1) - gap_extend
          _I_curr[row] = sub_sat(_I_prev[row], _gap_extend);  // I: I(i,j-1) - gap_extend
          // _tmp0 = S(i,j-1) - (gap_open + gap_extend)
          // Sc is _S_curr[col - 1], seed->S_col[row] for col=0
          _tmp0 = sub_sat(Sc, _gap_extend + _gap_open);
          _I_curr[row] = max(_I_curr[row], _tmp0);

          PROF_END("I")
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

          PROF_BEGIN("M")

          _Ceq = ZERO_CT;
          _Cneq = ZERO_CT;

          if (ref != Base::N) {
              // Set all mismatching pairs to _mismatch
              _tmp0 = cmp_neq(read, ref);
              _Cneq = _tmp0 & _mismatch;   // If the read base is Base::N, set to 0 (_Ceq)
              _tmp0 = cmp_eq(read, Base::N);
              _Cneq = blend(_Ceq, _Cneq, _tmp0);

              // b is not N, so all equal bases are valid
              _tmp0 = cmp_eq(read, ref);
              _Ceq = _tmp0 & _match;
          }

          // Sp is _S_prev[col - 1], 0 for row=0
          // Seed->S_col[row - 1] for col=0
          _S_curr[col] = add_sat(Sp, _Ceq);   // Add match scores
          _S_curr[col] = sub_sat(_S_curr[col], _Cneq); // Subtract mismatch scores

          PROF_END("M")
      }


      /**
       * @brief
       * Takes the max of D,I, and M vectors and stores the _curr_posent best score/position
       * Currently does not support non-default template args
       * @param row _curr_posent row
       * @param col _curr_posent column
       * @param node_origin Current position, used to get absolute alignment position
       */
      __INLINE__
      void _fill_cell_finish(const uint32_t &row,
                             const uint32_t &col,
                             const uint32_t &node_origin,
                             const std::vector<uint32_t> &targets) {
          using namespace simdpp;

          _curr_pos = node_origin + col + 1;    // absolute position in reference sequence

          // S(i,j) = max{ D(i,j), I(i,j), S(i-1,j-1) + C(s,t) }
          _S_curr[col] = max(_D_curr[col], _S_curr[col]);
          _S_curr[col] = max(_I_curr[row], _S_curr[col]);

          // Check for new or equal high scores
          PROF_BEGIN("New Max")

          _tmp0 = _S_curr[col] > _max_score;
          if (reduce_or(_tmp0)) {
              _max_score = max(_max_score, _S_curr[col]);

              for (uint8_t i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  // Check if i'th mask element is set
                  if (simdpp::extract<0>(_tmp0)) {
                      _max_elem = extract(i, _max_score);
                      // If current position is far enough away from the old max, demote the old max to sub_max
                      if (_curr_pos > _max_pos[i] + _read_len) {
                          insert(_max_elem, i, _sub_score);
                          _sub_pos[i] = _max_pos[i];
                          _sub_count[i] = _max_count[i];
                          if (_cor_flag[i] == 1) _cor_flag[i] = 2;
                      }
                      // max pos updated in eq check, same with cor flag
                      _max_count[i] = 0;
                  }
                  _tmp0 = simdpp::move16_l<1>(_tmp0);
              }
          }

          PROF_END("New Max")
          PROF_BEGIN("Eq Max")

          // Check for equal max score. If we set a new high score this will set the count to 1
          _tmp0 = cmp_eq(_S_curr[col], _max_score);
          if (reduce_or(_tmp0)) {
              for (uint8_t i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  // Check if the i'th elements MSB is set
                  if (simdpp::extract<0>(_tmp0)) {
                      _max_count[i] += 1;
                      _max_pos[i] = _curr_pos;
                      if (_curr_pos == targets[i]) _cor_flag[i] = 1;
                  }
                  _tmp0 = move16_l<1>(_tmp0);
              }
          }

          PROF_END("Eq Max")
          PROF_BEGIN("New Sub")

          // new second best score
          _tmp0 = _S_curr[col] > _sub_score;
          if (reduce_or(_tmp0)) {
              for (uint8_t i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  // Check if the i'th elements MSB is set, and if far enough
                  if (simdpp::extract<0>(_tmp0) && _max_pos[i] < _curr_pos - _read_len) {
                      //TODO add -(_max_score/_sub_score) term
                      insert(extract(i, _S_curr[col]), i, _sub_score);
                      // sub_pos will be set in eq check, same with cor flag
                      _sub_count[i] = 0;
                  }
                  _tmp0 = move16_l<1>(_tmp0);
              }
          }

          PROF_END("New Sub")
          PROF_BEGIN("Eq Sub")

          // Repeat sub score
          _tmp0 = cmp_eq(_S_curr[col], _sub_score);
          if (reduce_or(_tmp0)) {
              for (uint8_t i = 0; i < SIMDPP_FAST_INT8_SIZE; ++i) {
                  // Check if the i'th elements MSB is set
                  if (simdpp::extract<0>(_tmp0) && _max_pos[i] < _curr_pos - _read_len) {
                      _sub_count[i] += 1;
                      _sub_pos[i] = _curr_pos;
                      if (_curr_pos == targets[i] && !_cor_flag[i]) _cor_flag[i] = 2;
                  }
                  _tmp0 = move16_l<1>(_tmp0);
              }
          }

          PROF_END("Eq Sub")
      }

      /**
       * @brief
       * Extract the i'th element from a vector. No range checking is done.
       * @param i index of element
       * @param vec vector to extract from
       */
      __INLINE__
      uint8_t extract(uint8_t i, const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &vec) {
          return ((uint8_t *) &vec)[i];
      }

      /**
       * @brief
       * Insert into the i'th element from a vector. No range checking is done.
       * @param elem element to insert
       * @param i index of element
       * @param vec vector to insert in
       */
      __INLINE__
      void insert(uint8_t elem, uint8_t i, const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &vec) {
          ((uint8_t *) &vec)[i] = elem;
      }

      /*********************************** Variables ***********************************/

      AlignmentGroup _alignment_group;
      size_t _read_len; /**< Maximum read length. */

      // Zero vector
      const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> ZERO_CT = simdpp::splat(0);

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
          _tmp0,  /**< temporary for use within functions */
          _Ceq,  /**< Match score when read_base == ref_base */
          _Cneq;
      /**< mismatch penalty */

      uint8_t _max_elem;
      uint32_t _curr_pos;

      // Optimal alignment info
      simdpp::uint8<SIMDPP_FAST_INT8_SIZE> _max_score;
      uint32_t *_max_pos;
      uint8_t *_max_count;

      // Suboptimal alignment info
      simdpp::uint8<SIMDPP_FAST_INT8_SIZE> _sub_score;
      uint32_t *_sub_pos;
      uint8_t *_sub_count;

      uint8_t *_cor_flag;

      size_t _max_node_len;

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
