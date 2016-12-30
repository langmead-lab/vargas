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

#include "scoring.h"
#include "utils.h"
#include "simdpp/simd.h"
#include "graph.h"
#include "doctest.h"

#include <array>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <string>

#define VA_ALIGN_DEBUG 0 // Print score matrix if 1
#define SW_GRID 0 // SIMD element to print for debug

#if VA_ALIGN_DEBUG
#include <iomanip>
#endif

#define DEFAULT_TOL_FACTOR 4 // If the pos is +- read_len/tol, count as correct alignment

namespace vargas {

  /**
   * Corresponding from simdpp type from native type, where NATIVE_T = {uint8_t, uint16_t, uint32_t}
   */
  template <typename NATIVE_T>
  using SIMD_T = typename std::conditional<std::is_same<uint8_t, NATIVE_T>::value,
                                           simdpp::uint8<SIMDPP_FAST_INT8_SIZE>,
                                           typename std::conditional<std::is_same<uint16_t, NATIVE_T>::value,
                                                                     simdpp::uint16<SIMDPP_FAST_INT16_SIZE>,
                                                                     simdpp::uint32<SIMDPP_FAST_INT32_SIZE>>::type>::type;

  /**
* @brief
* Extract the i'th element from a vector. No range checking is done.
* @param i index of element
* @param vec vector to extract from
*/
  template <typename NATIVE_T>
  __RG_STRONG_INLINE__
  constexpr NATIVE_T extract(const uint8_t i, const SIMD_T<NATIVE_T> &vec) {
      return reinterpret_cast<const NATIVE_T *>(&vec)[i];
  }

  /**
   * @brief
   * Insert into the i'th element from a vector. No range checking is done.
   * @param elem element to insert
   * @param i index of element
   * @param vec vector to insert in
   */
  template <typename NATIVE_T>
  __RG_STRONG_INLINE__
  void insert(const NATIVE_T elem, const uint8_t i, SIMD_T<NATIVE_T> &vec) {
      reinterpret_cast<NATIVE_T *>(&vec)[i] = elem;
  }

  /**
   * @brief
   * Only makes sense for mask types, in base case MSB of each byte is checked.
   * @param c
   * @return Return true if c != 0
   */
  __RG_STRONG_INLINE__
  bool any_set(const simdpp::uint8x16 &c) {
      return simdpp::extract_bits_any(c);
  }

  __RG_STRONG_INLINE__
  bool any_set(const simdpp::uint16x8 &c) {
      return extract_bits_any(simdpp::bit_cast<simdpp::uint8x16, simdpp::uint16x8>(c));
  }

  __RG_STRONG_INLINE__
  bool any_set(const simdpp::uint32x4 &c) {
      return extract_bits_any(simdpp::bit_cast<simdpp::uint8x16, simdpp::uint32x4>(c));
  }

  __RG_STRONG_INLINE__
  bool any_set(const simdpp::uint8x32 &c) {
      return extract_bits_any(c);
  }
  __RG_STRONG_INLINE__
  bool any_set(const simdpp::uint16x16 &c) {
      return any_set(simdpp::bit_cast<simdpp::uint8x32, simdpp::uint16x16>(c));
  }

  __RG_STRONG_INLINE__
  bool any_set(const simdpp::uint16<32> &c) {
      simdpp::uint16x16 a, b;
      simdpp::split(c, a, b);
      return any_set(a) || any_set(b);
  }

  /**
   * @brief
   * Common base class to all template versions of Aligners
   */
  class AlignerBase {
    public:
      virtual ~AlignerBase() = 0;

      /**
       * Set the scoring scheme used for the alignments.
       * @param match match score
       * @param mismatch mismatch score
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      virtual void set_scores(uint8_t, uint8_t, uint8_t, uint8_t) = 0;

      /**
       * @brief
       * Set the score profile from a ScoreProfile - gap penalties between ref/read may vary.
       * @param prof
       */
      virtual void set_scores(const ScoreProfile &prof) = 0;

      /**
       * @brief
       * If the best score is +/- this tolerence of the target, count it as correct.
       * Impacts the corflag.
       * @param tol
       */
      virtual void set_correctness_tolerance(const size_t) = 0;

      /**
       * @return Current correctness tolerance
       */
      virtual size_t tolerance() const = 0;

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group vector of reads to align to
       * @param target high-score cell, determines correctness_flag, 1 based
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @param aligns Results packet to populate
       */
      virtual void align_into(const std::vector<std::string> &, const std::vector<size_t> &,
                              Graph::const_iterator, Graph::const_iterator, Results &) = 0;

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
      virtual Results align(const std::vector<std::string> &read_group, const std::vector<size_t> &targets,
                            Graph::const_iterator begin, Graph::const_iterator end) {
          Results aligns;
          align_into(read_group, targets, begin, end, aligns);
          return aligns;
      }

      /**
 * @brief
 * Align a batch of reads to a graph range, return a vector of alignments
 * corresponding to the reads.
 * @param read_group vector of reads to align to
 * @param begin iterator to beginning of graph
 * @param end iterator to end of graph
 * @return Results packet
 */
      virtual Results align(const std::vector<std::string> &read_group,
                            Graph::const_iterator begin, Graph::const_iterator end) {
          std::vector<size_t> targets(read_group.size());
          std::fill(targets.begin(), targets.end(), 0);
          return align(read_group, targets, begin, end);
      }

      static constexpr size_t default_tolerance() { return DEFAULT_TOL_FACTOR; }

    protected:
      ScoreProfile _prof;

  };
  inline AlignerBase::~AlignerBase() {}

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
   * Vargas::GraphFactory gb("reference.fa", "var.bcf");
   * gb.node_len(5);
   * gb.ingroup(100);
   * gb.region("x:0-15");
   *
   * Vargas::Graph g = gb.build();
   * std::vector<std::string> reads = {"ACGT", "GGGG", "ATTA", "CCNT"};
   *
   * Vargas::Aligner a(g.max_node_len(), 4);
   * Vargas::Aligner::Results res = a.align(reads, g.begin(), g.end());
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
   * @tparam NATIVE_T Native data type of score matrix element. One of uint8_t, uint16_t, uint32_t
   * @tparam END_TO_END If true, perform end to end alignment
   */
  template<typename NATIVE_T, bool END_TO_END = false>
  class AlignerT: public AlignerBase {
      static_assert(std::is_same<NATIVE_T, uint8_t>::value || std::is_same<NATIVE_T, uint16_t>::value, "Invalid type.");
    public:

      constexpr static size_t estimated_size(size_t nlen, size_t rlen) {
          return sizeof(AlignerT<NATIVE_T>)
          + (4 * sizeof(SIMD_T<NATIVE_T>) * nlen)
          + (2 * sizeof(SIMD_T<NATIVE_T>) * rlen);
      }

      /**
       * @brief
       * Set scoring parameters.
       * @param max_node_len max node length
       * @param read_len max read length
       * @param match match score
       * @param mismatch mismatch penalty
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      AlignerT(const size_t max_node_len, const size_t read_len,
               const uint8_t match = 2, const uint8_t mismatch = 2, const uint8_t open = 3, const uint8_t extend = 1) :
      _read_len(read_len),
      _max_node_len(max_node_len),
      _alignment_group(read_len) {
          set_scores(match, mismatch, open, extend);
          set_correctness_tolerance(read_len / DEFAULT_TOL_FACTOR);
          _alloc();
      }

      AlignerT(size_t max_node_len, size_t read_len, const ScoreProfile &prof) :
      _read_len(read_len),
      _max_node_len(max_node_len),
      _alignment_group(read_len) {
          set_scores(prof); // May throw
          set_correctness_tolerance(read_len / DEFAULT_TOL_FACTOR);
          _alloc();
      }


      ~AlignerT() noexcept {
          _dealloc();
      }

      AlignerT(const AlignerT<NATIVE_T, END_TO_END> &a) = delete;
      AlignerT(AlignerT<NATIVE_T, END_TO_END> &&a) = delete;
      AlignerT &operator=(const AlignerT<NATIVE_T> &) = delete;
      AlignerT &operator=(AlignerT<NATIVE_T> &&) = delete;

      /**
       * @brief
       * Container for a packaged batch of reads.
       * @details
       * Reads are interleaved so each SIMD vector
       * contains bases from all reads, respective to the base number. For example AlignmentGroup[0]
       * would contain the first bases of every read. All reads must be the same length. Minimal error checking.
       */
      class AlignmentGroup {
        public:

          AlignmentGroup(size_t read_len) : _read_len(read_len), _packaged_reads(read_len) {}

          __RG_STRONG_INLINE__
          void load_reads(const std::vector<std::string> &reads, const size_t begin, const size_t end) {
              load_reads(std::vector<std::string>(reads.begin() + begin, reads.begin() + end));
          }

          /**
           * @param batch load the given vector of reads.
           */
          __RG_STRONG_INLINE__
          void load_reads(const std::vector<std::string> &batch) {
              std::vector<std::vector<rg::Base>> _reads;
              for (auto &b : batch) _reads.push_back(rg::seq_to_num(b));
              load_reads(_reads);
          }

          /**
           * @param batch load the given vector of reads.
           */
          void load_reads(const std::vector<std::vector<rg::Base>> &batch) {
              _package_reads(batch);
          }

          /**
           * @brief
           * Return the i'th base of every read in a simdpp vector.
           * @param i base index.
           */
          const SIMD_T<NATIVE_T> &at(const size_t i) const {
              return _packaged_reads.at(i);
          }

          /**
           * @brief
           * Pointer to raw packaged read data.
           * @return VEC_TYPE pointer
           */
          const SIMD_T<NATIVE_T> *data() const {
              return _packaged_reads.data();
          }

          /**
           * @brief
           * Non const version of at(i).
           * @param i base index
           */
          SIMD_T<NATIVE_T> &operator[](const int i) {
              return _packaged_reads[i];
          }

          /**
           * @brief
           * Returns optimal number of reads in a batch based on SIMD architecture.
           * @return batch size.
           */
          static constexpr size_t group_size() { return SIMD_T<NATIVE_T>::length; }

          /**
           * @return iterator to the beginning of the packaged reads.
           */
          typename std::vector<SIMD_T<NATIVE_T>>::const_iterator begin() const {
              return _packaged_reads.cbegin();
          }

          /**
           * @return iterator to the end of the packaged reads.
           */
          typename std::vector<SIMD_T<NATIVE_T>>::const_iterator end() const {
              return _packaged_reads.cend();
          }

        private:

          const size_t _read_len;

          /**
           * _packaged_reads[i] contains all i'th bases.
           * The length of _packaged_reads is the length of the read,
           * where as the length of _packaged_reads[i] is the number
           * of reads.
           */
          std::vector<SIMD_T<NATIVE_T>> _packaged_reads;

          /**
           * Interleaves reads so all same-index base positions are in one
           * vector. Empty spaces are padded with Base::N.
           * @param _reads vector of reads to package
           */
          __RG_STRONG_INLINE__
          void _package_reads(const std::vector<std::vector<rg::Base>> &_reads) {
              assert(_reads.size() <= group_size());
              // Interleave reads
              // For each read (read[i] is in _packaged_reads[0..n][i]
              for (size_t r = 0; r < _reads.size(); ++r) {
                  assert(_reads[r].size() == _read_len);
                  // Put each base in the appropriate vector element
                  for (size_t p = 0; p < _read_len; ++p) {
                      insert<NATIVE_T>(_reads[r][p], r, _packaged_reads[p]);
                  }
              }

              // Pad underful batches
              for (size_t r = _reads.size(); r < group_size(); ++r) {
                  for (size_t p = 0; p < _read_len; ++p) {
                      insert<NATIVE_T>(rg::Base::N, r, _packaged_reads[p]);
                  }
              }
          }

      };

      void set_scores(uint8_t match, uint8_t mismatch, uint8_t open, uint8_t extend) override {
          _prof = ScoreProfile(match, mismatch, open, extend);
          _prof.end_to_end = END_TO_END;
          _match_vec = simdpp::splat<SIMD_T<NATIVE_T>>(match);
          _mismatch_vec = simdpp::splat<SIMD_T<NATIVE_T>>(mismatch);
          _gap_open_extend_vec_rd = simdpp::splat<SIMD_T<NATIVE_T>>(open + extend);
          _gap_extend_vec_rd = simdpp::splat<SIMD_T<NATIVE_T>>(extend);
          _gap_open_extend_vec_ref = simdpp::splat<SIMD_T<NATIVE_T>>(open + extend);
          _gap_extend_vec_ref = simdpp::splat<SIMD_T<NATIVE_T>>(extend);
          _ambig_vec = ZERO_CT;
          _bias = _get_bias(_read_len, match, mismatch, open, extend);
          _bias_vec = simdpp::splat<SIMD_T<NATIVE_T>>(_bias);
      }

      void set_scores(const ScoreProfile &prof) override {
          _prof = prof;
          _prof.end_to_end = END_TO_END;
          _match_vec = simdpp::splat<SIMD_T<NATIVE_T>>(prof.match);
          _mismatch_vec = simdpp::splat<SIMD_T<NATIVE_T>>(prof.mismatch);
          _gap_open_extend_vec_rd = simdpp::splat<SIMD_T<NATIVE_T>>(prof.read_gopen + prof.read_gext);
          _gap_extend_vec_rd = simdpp::splat<SIMD_T<NATIVE_T>>(prof.read_gext);
          _gap_open_extend_vec_ref = simdpp::splat<SIMD_T<NATIVE_T>>(prof.ref_gopen + prof.ref_gext);
          _gap_extend_vec_ref = simdpp::splat<SIMD_T<NATIVE_T>>(prof.ref_gext);
          _ambig_vec = simdpp::splat<SIMD_T<NATIVE_T>>(prof.ambig);
          _bias = _get_bias(_read_len, prof.match, prof.mismatch, prof.read_gopen, prof.read_gext);
          _bias_vec = simdpp::splat<SIMD_T<NATIVE_T>>(_bias);
      }

      void set_correctness_tolerance(const size_t tol) override {
          _prof.tol = tol;
      }

      size_t tolerance() const override {
          return _prof.tol;
      }

      /**
       * @return maximum number of reads that can be aligned at once.
       */
      static constexpr size_t read_capacity() { return SIMD_T<NATIVE_T>::length; }

      void align_into(const std::vector<std::string> &read_group, const std::vector<size_t> &targets,
                      Graph::const_iterator begin, Graph::const_iterator end, Results &aligns) override {

          assert(targets.size() == read_group.size());

          const size_t num_groups = 1 + ((read_group.size() - 1) / read_capacity());
          // Possible oversize if there is a partial group
          aligns.resize(num_groups * read_capacity());
          _targets_lower.resize(num_groups * read_capacity());
          _targets_upper.resize(num_groups * read_capacity());

          std::transform(targets.begin(), targets.end(), _targets_lower.begin(),
                         [this](size_t x) { return x - _prof.tol; });
          std::transform(targets.begin(), targets.end(), _targets_upper.begin(),
                         [this](size_t x) { return x + _prof.tol; });
          std::fill(aligns.correct.begin(), aligns.correct.end(), 0);

          // Keep the scores at the positions, overwrites position. [0] is current position, 1-:ead_capacity + 1 is pos
          std::unordered_map<size_t, _seed> seed_map; // Maps node ID's to the ending matrix columns of the node
          _seed seed(_read_len);

          for (size_t group = 0; group < num_groups; ++group) {
              seed_map.clear();

              // Subset of read set
              const size_t beg_offset = group * read_capacity();
              const size_t end_offset = std::min((group + 1) * read_capacity(), read_group.size());
              const size_t len = end_offset - beg_offset;
              assert(len <= read_capacity());

              _alignment_group.load_reads(read_group, beg_offset, end_offset);
              _max_score = simdpp::splat<SIMD_T<NATIVE_T>>(std::numeric_limits<NATIVE_T>::min());
              _sub_score = _max_score;
              _max_pos = aligns.max_pos.data() + beg_offset;
              _sub_pos = aligns.sub_pos.data() + beg_offset;
              _max_count = aligns.max_count.data() + beg_offset;
              _sub_count = aligns.sub_count.data() + beg_offset;
              _cor_flag = aligns.correct.data() + beg_offset;
              _targets_lower_ptr = _targets_lower.data() + beg_offset;
              _targets_upper_ptr = _targets_upper.data() + beg_offset;

              // Create subrange of targets, and sort by position.
              for (uint8_t j = 0; j < len; ++j) {
                  _target_subrange[j].idx = j;
                  _target_subrange[j].pos = targets[beg_offset + j];
                  _target_subrange[j].score = std::numeric_limits<int>::min();
              }
              for (size_t j = len; j < read_capacity(); ++j) { _target_subrange[j].idx = -1; }
              std::sort(_target_subrange.begin(), _target_subrange.end(),
                        [](const _target &p, const _target &q) { return q.idx < 0 || (p.idx >= 0 && p.pos < q.pos); });

              // seed matrix and first node
              _seed_matrix(seed);
              _fill_node(*begin, _alignment_group, seed, seed_map.emplace(begin->id(), _read_len).first->second);
              for (auto gi = begin + 1; gi != end; ++gi) {
                  _get_seed(gi.incoming(), seed_map, seed);
                  if (gi->is_pinched()) seed_map.clear();
                  _fill_node(*gi, _alignment_group, seed, seed_map.emplace(gi->id(), _read_len).first->second);
              }

              // Copy scores
              for (uint8_t i = 0; i < len; ++i) {
                  aligns.max_score[beg_offset + i] = int(vargas::extract<NATIVE_T>(i, _max_score)) - _bias;
                  aligns.sub_score[beg_offset + i] = int(vargas::extract<NATIVE_T>(i, _sub_score)) - _bias;
                  aligns.target_score[beg_offset + _target_subrange[i].idx] = int(_target_subrange[i].score) - _bias;

              }
          }
          // Crop off potential buffer
          aligns.resize(read_group.size());
          aligns.profile = _prof;
      }

    private:

      /**
       * @brief
       * Ending vectors from a previous node
       */
      struct _seed {
          _seed() = delete;
          _seed(const size_t _read_len) : S_col(_read_len) {}
          std::vector<SIMD_T<NATIVE_T>> S_col; /**< Last column of score matrix.*/
      };

      struct _target {
          int idx;
          int score;
          size_t pos;
      };

      void _seed_matrix(_seed &seed) {
          std::fill(seed.S_col.begin(), seed.S_col.end(), _bias_vec);

          if (END_TO_END) {
              for (size_t i = 0; i < _read_len; ++i) {
                  const uint8_t b = _prof.read_gopen + (i * _prof.read_gext);
                  seed.S_col[i] = simdpp::sub_sat(seed.S_col[i], b);
              }
          }
      }

      /**
       * @brief
       * Returns the best seed from all previous nodes.
       * Graph should be validated before alignment to ensure proper seed fetch
       * @param prev_ids All nodes preceding _curr_posent node
       * @param seed_map ID->seed map for all previous nodes
       * @param seed best seed to populate
       * @throws std::domain_error if a node listed as a previous node but it has not been encountered yet.
       * i.e. not topographically sorted.
       */
      __RG_STRONG_INLINE__
      void _get_seed(const std::vector<size_t> &prev_ids,
                     std::unordered_map<size_t, _seed> &seed_map, _seed &seed) const {
          try {
              for (size_t i = 0; i < _read_len; ++i) {
                  seed.S_col[i] = ZERO_CT;
                  for (size_t id : prev_ids) {
                      seed.S_col[i] = simdpp::max(seed.S_col[i], seed_map.at(id).S_col[i]);
                  }
              }
          }
          catch (std::exception &e) {
              throw std::domain_error("Invalid node ordering.");
          }
      }

      __RG_STRONG_INLINE__
      void _set_nullptr() noexcept {
          _Sa = _Sb = _Da = _Db = _Ia = _Ib = nullptr;
          _S_prev = _S_curr = _D_prev = _D_curr = _I_prev = _I_curr = nullptr;
      }

      /**
       * @brief
       * deletes allocated matrix filling vectors.
       */
      void _dealloc() noexcept {
          if (_Sa) delete _Sa;
          if (_Sb) delete _Sb;
          if (_Da) delete _Da;
          if (_Db) delete _Db;
          if (_Ia) delete _Ia;
          if (_Ib) delete _Ib;
          _set_nullptr();
      }

      /**
       * @brief
       * Allocate S and D vectors.
       */
      void _alloc() {
          _Sa = new std::vector<SIMD_T<NATIVE_T>>(_max_node_len);
          _Sb = new std::vector<SIMD_T<NATIVE_T>>(_max_node_len);
          _Da = new std::vector<SIMD_T<NATIVE_T>>(_max_node_len);
          _Db = new std::vector<SIMD_T<NATIVE_T>>(_max_node_len);
          _Ia = new std::vector<SIMD_T<NATIVE_T>>(_read_len);
          _Ib = new std::vector<SIMD_T<NATIVE_T>>(_read_len);

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
      __RG_STRONG_INLINE__
      void _fill_node(const Graph::Node &n, const AlignmentGroup &read_group, const _seed &s, _seed &nxt) {

          #if VA_ALIGN_DEBUG
          const size_t fw = _bias > 0 ? 4 : 3;
          std::cerr << "X " << std::setw(fw) << (int)_bias << "|";
          for(const char c : n.seq_str()) std::cerr << std::setw(fw) << c << " ";
          std::cerr << '\n';
          #endif

          // Empty nodes represents deletions
          if (n.seq().size() == 0) {
              nxt = s;
              return;
          }

          assert(n.seq().size() <= _max_node_len);

          const rg::Base *node_seq = n.seq().data();
          const SIMD_T<NATIVE_T> *read_ptr = read_group.data();
          const size_t seq_size = n.seq().size();
          const size_t node_origin = n.end_pos() - seq_size + 2;

          int csp_start = 0;
          while (_target_subrange[csp_start].pos < node_origin) ++csp_start;
          int csp = csp_start;

          // top left corner
          _fill_cell_rzcz(read_ptr[0], node_seq[0], s);
          _fill_cell_finish_prox(0, 0, node_origin, csp);

          // top row
          for (size_t c = 1; c < seq_size; ++c) {
              _fill_cell_rz(read_ptr[0], node_seq[c], c);
              _fill_cell_finish_prox(0, c, node_origin, csp);
          }

          #if VA_ALIGN_DEBUG
          std::cerr << rg::num_to_base((rg::Base)simdpp::extract<SW_GRID>(read_ptr[0])) << " " << std::setw(fw)
                    << (int)simdpp::extract<SW_GRID>(s.S_col[0]) << "|";
          for (size_t c = 0; c < seq_size; ++c) {
              std::cerr << std::setw(fw) << (int)simdpp::extract<SW_GRID>(_S_curr[c]) - _bias << " ";
          }
          std::cerr << "\n";
          #endif

          nxt.S_col[0] = _S_curr[seq_size - 1];

          // Rest of the rows
          for (size_t r = 1; r < _read_len; ++r) {
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

              csp = csp_start;

              // first col
              _fill_cell_cz(read_ptr[r], node_seq[0], r, s);
              _fill_cell_finish_prox(r, 0, node_origin, csp);

              // Inner grid
              for (size_t c = 1; c < seq_size; ++c) {
                  _fill_cell(read_ptr[r], node_seq[c], r, c);
                  _fill_cell_finish_prox(r, c, node_origin, csp);
              }

              #if VA_ALIGN_DEBUG
              std::cerr << rg::num_to_base((rg::Base)simdpp::extract<SW_GRID>(read_ptr[r])) << " " << std::setw(fw)
                        << (int)simdpp::extract<SW_GRID>(s.S_col[r]) << "|";
              for (size_t c = 0; c < seq_size; ++c) {
                  std::cerr << std::setw(fw) << (int)simdpp::extract<SW_GRID>(_S_curr[c]) - _bias << " ";
              }
              std::cerr << "\n";
              #endif

              nxt.S_col[r] = _S_curr[seq_size - 1];

          }

          #if VA_ALIGN_DEBUG
          std::cerr << "\n";
          #endif

          // Look at last row for best scores when end to end
          if (END_TO_END) {
              for (size_t c = 0; c < seq_size; ++c) {
                  _fill_cell_finish(c, node_origin, csp);
              }
          }

      }

      /**
       * @brief
       * Fills the top left cell.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param s alignment seed from previous node
       */
      __RG_STRONG_INLINE__
      void _fill_cell_rzcz(const SIMD_T<NATIVE_T> &read_base, const rg::Base &ref, const _seed &s) {
          _D(0, _bias_vec, _bias_vec);
          _I(0, s.S_col[0]);
          _M(0, read_base, ref, _bias_vec);
      }

      /**
       * @brief
       * Fills cells when row is 0.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param col _curr_posent column in matrix
       */
      __RG_STRONG_INLINE__
      void _fill_cell_rz(const SIMD_T<NATIVE_T> &read_base, const rg::Base &ref, const size_t &col) {
          _D(col, _bias_vec, _bias_vec);
          _I(0, _S_curr[col - 1]);
          _M(col, read_base, ref, _bias_vec);
      }

      /**
       * @brief
       * Fills cells when col is 0.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row _curr_posent row in matrix
       * @param s alignment seed from previous node
       */
      __RG_STRONG_INLINE__
      void _fill_cell_cz(const SIMD_T<NATIVE_T> &read_base, const rg::Base &ref, const size_t &row, const _seed &s) {
          _D(0, _D_prev[0], _S_prev[0]);
          _I(row, s.S_col[row]);
          _M(0, read_base, ref, s.S_col[row - 1]);
      }

      /**
       * @brief
       * Fills the _curr_posent cell.
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row _curr_posent row in matrix
       * @param col _curr_posent column in matrix
       */
      __RG_STRONG_INLINE__
      void _fill_cell(const SIMD_T<NATIVE_T> &read_base, const rg::Base &ref, const size_t &row, const size_t &col) {
          _D(col, _D_prev[col], _S_prev[col]);
          _I(row, _S_curr[col - 1]);
          _M(col, read_base, ref, _S_prev[col - 1]);
      }

      /**
       * @brief
       * Score if there is a deletion in read (gap in read)
       * @param col _curr_posent column
       * @param Dp Previous D value at _curr_posent col.
       * @param Sp Previous S value at _curr_posent col.
       */
      __RG_STRONG_INLINE__
      void _D(const size_t &col, const SIMD_T<NATIVE_T> &Dp, const SIMD_T<NATIVE_T> &Sp) {
          using namespace simdpp;

          // D(i,j) = D(i-1,j) - gap_extend
          // Dp is _D_prev[col], 0 for row=0
          _D_curr[col] = sub_sat(Dp, _gap_extend_vec_rd);   // _tmp0 = S(i-1,j) - ( gap_open + gap_extend)
          // Sp is _S_prev[col], 0 for row=0
          // D(i,j) = max{ D(i-1,j) - gap_extend, S(i-1,j) - ( gap_open + gap_extend) }
          _tmp0 = sub_sat(Sp, _gap_open_extend_vec_rd);
          _D_curr[col] = max(_D_curr[col], _tmp0);

      }

      /**
       * @brief
       * Score if there is an insertion in read (gap in ref)
       * @param row _curr_posent row
       * @param Sc Previous S value (cell to the left)
       */
      __RG_STRONG_INLINE__
      void _I(const size_t &row, const SIMD_T<NATIVE_T> &Sc) {
          using namespace simdpp;

          // I(i,j) = I(i,j-1) - gap_extend
          _I_curr[row] = sub_sat(_I_prev[row], _gap_extend_vec_ref);  // I: I(i,j-1) - gap_extend
          // _tmp0 = S(i,j-1) - (gap_open + gap_extend)
          // Sc is _S_curr[col - 1], seed->S_col[row] for col=0
          _tmp0 = sub_sat(Sc, _gap_open_extend_vec_ref);
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
      __RG_STRONG_INLINE__
      void _M(size_t col, const SIMD_T<NATIVE_T> &read, const rg::Base &ref, const SIMD_T<NATIVE_T> &Sp) {
          using namespace simdpp;

          if (ref != rg::Base::N) {
              _Ceq = ZERO_CT;
              _Cneq = ZERO_CT;

              // Set all mismatching pairs to _mismatch
              _tmp0 = read != _base_vec[ref];
              _Cneq = _tmp0 & _mismatch_vec;
              // If the read base is Base::N, set to ambig penalty
              _tmp0 = read == _base_vec[rg::Base::N];
              _Cneq = blend(_ambig_vec, _Cneq, _tmp0);

              // b is not N, so all equal bases are valid
              _tmp0 = read == _base_vec[ref];
              _Ceq = _tmp0 & _match_vec;

              // Sp is _S_prev[col - 1], 0 for row=0
              // Seed->S_col[row - 1] for col=0
              _S_curr[col] = add_sat(Sp, _Ceq);   // Add match scores
              _S_curr[col] = sub_sat(_S_curr[col], _Cneq); // Subtract mismatch scores
          } else {
              // Penalty for N
              _S_curr[col] = sub_sat(Sp, _ambig_vec);
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
      __RG_STRONG_INLINE__ __RG_UNROLL__
      void _fill_cell_finish(const size_t &col, const size_t &node_origin, int &curr_search) {
          using namespace simdpp;

          _curr_pos = node_origin + col;    // absolute position in reference sequence
          while (_target_subrange[curr_search].pos == _curr_pos) {
              _target_subrange[curr_search].score =
              std::max<int>(vargas::extract<NATIVE_T>(_target_subrange[curr_search].idx, _S_curr[col]),
                            _target_subrange[curr_search].score);
              ++curr_search;
          }


          _tmp0 = _S_curr[col] > _max_score;
          if (any_set(_tmp0)) {
              // Check for new or equal high scores
              _max_score = max(_S_curr[col], _max_score);
              for (size_t i = 0; i < read_capacity(); ++i) {
                  if (_tmp0_ptr[i]) {
                      // Demote old max to submax
                      if (_curr_pos > _max_pos[i] + _read_len) {
                          _sub_score_ptr[i] = _max_score_ptr[i];
                          _sub_pos[i] = _max_pos[i];
                          _sub_count[i] = _max_count[i];
                          if (_cor_flag[i] == 1) _cor_flag[i] = 2;
                          else _cor_flag[i] = 0;
                      }
                      _max_count[i] = 1;
                      if (_cor_flag[i] == 1) _cor_flag[i] = 0;
                      else _max_pos[i] = _curr_pos;
                  }
              }
          }

          _tmp0 = cmp_eq(_S_curr[col], _max_score);
          if (any_set(_tmp0)) {
              // Check for equal max score.
              for (size_t i = 0; i < read_capacity(); ++i) {
                  if (_tmp0_ptr[i]) {
                      if (_curr_pos > _max_pos[i] + _read_len) ++(_max_count[i]);
                      _max_pos[i] = _curr_pos;
                      if (_curr_pos >= _targets_lower_ptr[i] && _curr_pos <= _targets_upper_ptr[i]) _cor_flag[i] = 1;
                  }
              }
          }


          _tmp0 = cmp_eq(_S_curr[col], _sub_score);
          if (any_set(_tmp0)) {
              // Repeat sub score
              for (size_t i = 0; i < read_capacity(); ++i) {
                  if (_tmp0_ptr[i] && _curr_pos > _max_pos[i] + _read_len) {
                      _sub_count[i] += _curr_pos > (_sub_pos[i] + _read_len);
                      _sub_pos[i] = _curr_pos;
                      if (_curr_pos >= _targets_lower_ptr[i] && _curr_pos <= _targets_upper_ptr[i]) _cor_flag[i] = 2;
                  }
              }
          }

          // Greater than old sub max and less than max score (prevent repeats of max triggering)
          _tmp0 = (_S_curr[col] > _sub_score) & (_S_curr[col] < _max_score);
          if (any_set(_tmp0)) {
              // new second best score
              for (size_t i = 0; i < read_capacity(); ++i) {
                  if (_tmp0_ptr[i] && _curr_pos > _max_pos[i] + _read_len) {
                      _sub_score_ptr[i] = vargas::extract<NATIVE_T>(i, _S_curr[col]);
                      _sub_count[i] = 1;
                      _sub_pos[i] = _curr_pos;
                      if (_curr_pos >= _targets_lower_ptr[i] && _curr_pos <= _targets_upper_ptr[i]) _cor_flag[i] = 2;
                      else _cor_flag[i] = _cor_flag[i] == 1;
                  }
              }
          }

      }

      /**
       * @brief
       * Chooses score selection based on end-to-end
       * @param row _curr_posent row
       * @param col _curr_posent column
       * @param node_origin Current position, used to get absolute alignment position
       */
      __RG_STRONG_INLINE__
      void _fill_cell_finish_prox(const size_t &row, const size_t &col, const size_t &node_origin,
                                  int &curr_search) {
          // S(i,j) = max{ D(i,j), I(i,j), S(i-1,j-1) + C(s,t) }
          _S_curr[col] = max(_D_curr[col], _S_curr[col]);
          _S_curr[col] = max(_I_curr[row], _S_curr[col]);

          if (!END_TO_END) _fill_cell_finish(col, node_origin, curr_search);
      }

      /**
       * @brief
       * Map each base to a vector of bases. Prevents repeated splat()
       */
      static const std::array<SIMD_T<NATIVE_T>, 5> _make_base_vec() {
          using rg::Base;
          static_assert(Base::A < 5 && Base::C < 5 && Base::G < 5  && Base::T < 5 && Base::N < 5, "Base enum error.");
          std::array<SIMD_T<NATIVE_T>, 5> v;
          v[Base::A] = simdpp::splat<SIMD_T<NATIVE_T>>(Base::A);
          v[Base::C] = simdpp::splat<SIMD_T<NATIVE_T>>(Base::C);
          v[Base::G] = simdpp::splat<SIMD_T<NATIVE_T>>(Base::G);
          v[Base::T] = simdpp::splat<SIMD_T<NATIVE_T>>(Base::T);
          v[Base::N] = simdpp::splat<SIMD_T<NATIVE_T>>(Base::N);
          return v;
      }

      static NATIVE_T _get_bias(const size_t read_len, const uint8_t match, const uint8_t mismatch,
                                const uint8_t gopen, const uint8_t gext) {
          static bool has_warned = false;
          if (read_len * match > std::numeric_limits<NATIVE_T>::max()) {
              throw std::domain_error("Insufficient bit-width for given match score and read length.");
          }
          if (!END_TO_END) return 0;

          // End to end alignment
          unsigned int b = std::numeric_limits<NATIVE_T>::max() - (read_len * match);

          //TODO Could be relaxed - all indels or all mismatch is unreasonable
          if (!has_warned && (gopen + (gext * (read_len - 1)) > b || read_len * mismatch > b)) {
              std::cerr << "WARN: Possibility of score saturation with parameters in end-to-end mode:\n\t"
                        << "Cell Width: " << (int) std::numeric_limits<NATIVE_T>::max() << ", "
                        << "Bias: " << b << "\n";
              has_warned = true;
          }
          return b;
      }


      /*********************************** Variables ***********************************/

      const size_t _read_len; /**< Maximum read length. */

      const std::array<SIMD_T<NATIVE_T>, 5> _base_vec = _make_base_vec();

      const SIMD_T<NATIVE_T> ZERO_CT = simdpp::splat<SIMD_T<NATIVE_T>>(0);
      NATIVE_T _bias;
      SIMD_T<NATIVE_T> _bias_vec;

      SIMD_T<NATIVE_T> _match_vec, _mismatch_vec,
      _gap_open_extend_vec_rd, _gap_extend_vec_rd, _gap_open_extend_vec_ref, _gap_extend_vec_ref, _ambig_vec;

      SIMD_T<NATIVE_T>
      *_S_prev = nullptr,    /**< _S_prev[n] => S(i-1, n) */
      *_S_curr = nullptr,
      *_D_prev = nullptr,    /**< _D_prev[n] => D(i-1, n) */
      *_D_curr = nullptr,
      *_I_prev = nullptr,    /**< _I_prev[r] => I(r, j-1) */
      *_I_curr = nullptr,
      *_swp_tmp0;

      /**
       * Each vector has an 'a' and a 'b' version. Through each row of the
       * matrix fill, their roles are swapped such that one becomes the previous
       * loops data, and the other is filled in.
       * S and D are padded 1 to provide a left column buffer.
       */
      std::vector<SIMD_T<NATIVE_T>>
      *_Sa = nullptr,    /**< Matrix row */
      *_Sb = nullptr,
      *_Da = nullptr,    /**< Deletion vector */
      *_Db = nullptr,
      *_Ia = nullptr,    /**< Insertion vector */
      *_Ib = nullptr;

      SIMD_T<NATIVE_T>
      _tmp0, /**< temporary for use within functions */
      _Ceq,  /**< Match score when read_base == ref_base */
      _Cneq, /**< mismatch penalty */
      _N_pen; /**< N penalty */

      size_t _curr_pos;

      // alignment info
      SIMD_T<NATIVE_T> _max_score;
      SIMD_T<NATIVE_T> _sub_score;
      size_t *_max_pos;
      size_t *_max_count;
      size_t *_sub_pos;
      size_t *_sub_count;

      NATIVE_T *const _tmp0_ptr = (NATIVE_T *) &_tmp0;
      NATIVE_T *const _max_score_ptr = (NATIVE_T *) &_max_score;
      NATIVE_T *const _sub_score_ptr = (NATIVE_T *) &_sub_score;

      uint8_t *_cor_flag;
      std::vector<size_t> _targets_lower, _targets_upper;
      size_t *_targets_lower_ptr, *_targets_upper_ptr;
      std::array<_target, SIMD_T<NATIVE_T>::length> _target_subrange;

      const size_t _max_node_len;

      AlignmentGroup _alignment_group;

  };

  using Aligner = AlignerT<uint8_t, false>;
  using WordAligner = AlignerT<uint16_t, false>;
  using AlignerETE = AlignerT<uint8_t, true>;
  using WordAlignerETE = AlignerT<uint16_t, true>;

}

TEST_SUITE("Aligners");

TEST_CASE ("Alignment") {

    vargas::Graph::Node::_newID = 0;
    vargas::Graph g;

    /**
    *     GGG
    *    /   \
    * AAA     TTTA
    *    \   /
    *     CCC(ref)
    */

    {
        vargas::Graph::Node n;
        n.set_endpos(2);
        n.set_as_ref();
        std::vector<bool> a = {0, 1, 1};
        n.set_population(a);
        n.set_seq("AAA");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(5);
        n.set_as_ref();
        std::vector<bool> a = {0, 0, 1};
        n.set_population(a);
        n.set_af(0.4);
        n.set_seq("CCC");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(5);
        n.set_not_ref();
        std::vector<bool> a = {0, 1, 0};
        n.set_population(a);
        n.set_af(0.6);
        n.set_seq("GGG");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
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

    SUBCASE("Graph Alignment") {
        std::vector<std::string> reads;
        reads.push_back("NNNCCTT");
        reads.push_back("NNNGGTT");
        reads.push_back("NNNAAGG");
        reads.push_back("NNNAACC");
        reads.push_back("NNAGGGT");
        reads.push_back("NNNNNGG");
        reads.push_back("AAATTTA");
        reads.push_back("AAAGCCC");
        const std::vector<size_t> origins = {8, 8, 5, 5, 7, 6, 10, 6};

        vargas::Aligner a(5, 7);
        vargas::Results aligns = a.align(reads, origins, g.begin(), g.end());
        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);
        CHECK((int) aligns.correct[0] == 1);
        CHECK(aligns.max_score[0] == aligns.target_score[0]); // since cf = 1

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);
        CHECK((int) aligns.correct[1] == 1);
        CHECK(aligns.max_score[1] == aligns.target_score[1]); // since cf = 1

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);
        CHECK((int) aligns.correct[2] == 1);
        CHECK(aligns.max_score[2] == aligns.target_score[2]); // since cf = 1

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);
        CHECK((int) aligns.correct[3] == 1);
        CHECK(aligns.max_score[3] == aligns.target_score[3]); // since cf = 1

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);
        CHECK((int) aligns.correct[4] == 1);
        CHECK(aligns.max_score[4] == aligns.target_score[4]); // since cf = 1

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);
        CHECK((int) aligns.correct[5] == 1);
        CHECK(aligns.max_score[5] == aligns.target_score[5]); // since cf = 1

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);
        CHECK((int) aligns.correct[6] == 1);
        CHECK(aligns.max_score[6] == aligns.target_score[6]); // since cf = 1

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);
        CHECK((int) aligns.correct[7] == 1);
        CHECK(aligns.max_score[7] == aligns.target_score[7]); // since cf = 1
    }

    SUBCASE("Scoring Scheme") {

        std::vector<std::string> reads;
        reads.push_back("NNNNNNCCTT");
        reads.push_back("NNNNNNGGTT");
        reads.push_back("NNNNNNAAGG");
        reads.push_back("NNNNNNAACC");
        reads.push_back("NNNNNAGGGT");
        reads.push_back("NNNNNNNNGG");
        reads.push_back("NNNAAATTTA");
        reads.push_back("NNNAAAGCCC");
        reads.push_back("AAAGAGTTTA");
        reads.push_back("AAAGAATTTA");
        const std::vector<size_t> origins = {8, 8, 5, 5, 7, 6, 10, 4, 10, 10};

        // hisat like params
        vargas::Aligner a(5, 10, 2, 6, 5, 3);
        vargas::Results aligns = a.align(reads, origins, g.begin(), g.end());

        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);
        CHECK((int) aligns.correct[0] == 1);
        CHECK(aligns.max_score[0] == aligns.target_score[0]);

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);
        CHECK((int) aligns.correct[1] == 1);
        CHECK(aligns.max_score[1] == aligns.target_score[1]);

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);
        CHECK((int) aligns.correct[2] == 1);
        CHECK(aligns.max_score[2] == aligns.target_score[2]);

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);
        CHECK((int) aligns.correct[3] == 1);
        CHECK(aligns.max_score[3] == aligns.target_score[3]);

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);
        CHECK((int) aligns.correct[4] == 1);
        CHECK(aligns.max_score[4] == aligns.target_score[4]);

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);
        CHECK((int) aligns.correct[5] == 1);
        CHECK(aligns.max_score[5] == aligns.target_score[5]);

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);
        CHECK((int) aligns.correct[6] == 1);
        CHECK(aligns.max_score[6] == aligns.target_score[6]);

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);
        CHECK((int) aligns.correct[7] == 1);
        CHECK(aligns.max_score[7] == aligns.target_score[7]);

        CHECK(aligns.max_score[8] == 12);
        CHECK(aligns.max_pos[8] == 10);
        CHECK((int) aligns.correct[8] == 1);
        CHECK(aligns.max_score[8] == aligns.target_score[8]);

        CHECK(aligns.max_score[9] == 8);
        CHECK(aligns.max_pos[9] == 10);
        CHECK((int) aligns.correct[9] == 1);
        CHECK(aligns.max_score[9] == aligns.target_score[9]);
    }

    SUBCASE("Scoring Scheme- N penalty") {

        std::vector<std::string> reads;
        reads.push_back("AAANGGTTTA");
        reads.push_back("AANNGGTTTA");
        reads.push_back("AAANNNTTTA");

        vargas::ScoreProfile prof(2, 2, 3, 1);
        prof.ambig = 1;
        vargas::Aligner a(5, 10, prof);
        vargas::Results aligns = a.align(reads, g.begin(), g.end());
        CHECK(aligns.max_score[0] == 17);
        CHECK(aligns.max_pos[0] == 10);

        CHECK(aligns.max_score[1] == 14);
        CHECK(aligns.max_pos[1] == 10);

        CHECK(aligns.max_score[2] == 11);
        CHECK(aligns.max_pos[2] == 10);
    }

    SUBCASE("Graph Alignment- Word") {
        std::vector<std::string> reads;
        reads.push_back("NNNCCTT");
        reads.push_back("NNNGGTT");
        reads.push_back("NNNAAGG");
        reads.push_back("NNNAACC");
        reads.push_back("NNAGGGT");
        reads.push_back("NNNNNGG");
        reads.push_back("AAATTTA");
        reads.push_back("AAAGCCC");
        const std::vector<size_t> origins = {8, 8, 5, 5, 7, 6, 10, 6};

        vargas::WordAligner a(5, 7);
        vargas::Results aligns = a.align(reads, origins, g.begin(), g.end());
        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);
        CHECK((int) aligns.correct[0] == 1);
        CHECK(aligns.max_score[0] == aligns.target_score[0]);

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);
        CHECK((int) aligns.correct[1] == 1);
        CHECK(aligns.max_score[1] == aligns.target_score[1]);

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);
        CHECK((int) aligns.correct[2] == 1);
        CHECK(aligns.max_score[2] == aligns.target_score[2]);

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);
        CHECK((int) aligns.correct[3] == 1);
        CHECK(aligns.max_score[3] == aligns.target_score[3]);

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);
        CHECK((int) aligns.correct[4] == 1);
        CHECK(aligns.max_score[4] == aligns.target_score[4]);

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);
        CHECK((int) aligns.correct[5] == 1);
        CHECK(aligns.max_score[5] == aligns.target_score[5]);

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);
        CHECK((int) aligns.correct[6] == 1);
        CHECK(aligns.max_score[6] == aligns.target_score[6]);

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);
        CHECK((int) aligns.correct[7] == 1);
        CHECK(aligns.max_score[7] == aligns.target_score[7]);
    }

    SUBCASE("Scoring Scheme- Word") {

        std::vector<std::string> reads;
        reads.push_back("NNNNNNCCTT");
        reads.push_back("NNNNNNGGTT");
        reads.push_back("NNNNNNAAGG");
        reads.push_back("NNNNNNAACC");
        reads.push_back("NNNNNAGGGT");
        reads.push_back("NNNNNNNNGG");
        reads.push_back("NNNAAATTTA");
        reads.push_back("NNNAAAGCCC");
        reads.push_back("AAAGAGTTTA");
        reads.push_back("AAAGAATTTA");
        const std::vector<size_t> origins = {8, 8, 5, 5, 7, 6, 10, 4, 10, 10};

        // hisat like params
        vargas::WordAligner a(5, 10, 2, 6, 5, 3);
        vargas::Results aligns = a.align(reads, origins, g.begin(), g.end());

        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);
        CHECK((int) aligns.correct[0] == 1);
        CHECK(aligns.max_score[0] == aligns.target_score[0]);

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);
        CHECK((int) aligns.correct[1] == 1);
        CHECK(aligns.max_score[1] == aligns.target_score[1]);

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);
        CHECK((int) aligns.correct[2] == 1);
        CHECK(aligns.max_score[2] == aligns.target_score[2]);

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);
        CHECK((int) aligns.correct[3] == 1);
        CHECK(aligns.max_score[3] == aligns.target_score[3]);

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);
        CHECK((int) aligns.correct[4] == 1);
        CHECK(aligns.max_score[4] == aligns.target_score[4]);

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);
        CHECK((int) aligns.correct[5] == 1);
        CHECK(aligns.max_score[5] == aligns.target_score[5]);

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);
        CHECK((int) aligns.correct[6] == 1);
        CHECK(aligns.max_score[6] == aligns.target_score[6]);

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);
        CHECK((int) aligns.correct[7] == 1);
        CHECK(aligns.max_score[7] == aligns.target_score[7]);

        CHECK(aligns.max_score[8] == 12);
        CHECK(aligns.max_pos[8] == 10);
        CHECK((int) aligns.correct[8] == 1);
        CHECK(aligns.max_score[8] == aligns.target_score[8]);

        CHECK(aligns.max_score[9] == 8);
        CHECK(aligns.max_pos[9] == 10);
        CHECK((int) aligns.correct[9] == 1);
        CHECK(aligns.max_score[9] == aligns.target_score[9]);
    }
}

TEST_CASE ("Indels") {
    vargas::Graph::Node::_newID = 0;
    vargas::Graph g;

    /**
     *ACTGCTNCAGTCAGTGNANACNCAC--ACGATCGTACGCNAGCTAGCCACAGTGCCCCCCTATATACGAN
     */

    {
        vargas::Graph::Node n;
        n.set_endpos(24);
        n.set_as_ref();
        n.set_seq("ACTGCTNCAGTCAGTGNANACNCAC");
        g.add_node(n);
    }

    {
        vargas::Graph::Node n;
        n.set_endpos(42);
        n.set_as_ref();
        n.set_seq("ACGATCGTACGCNAGCTAGCCACAGTGCCCCCCTATATACGAN");
        g.add_node(n);
    }
    g.add_edge(0, 1);

    SUBCASE ("Indel") {
        std::vector<std::string> reads;
        reads.push_back("ACTGCTNCAGTC"); // perfect alignment, pos 1
        reads.push_back("ACTGCTACAGTC"); // perfect alignment, pos 1, diff N
        reads.push_back("ACGATCACGCNA"); // 2 del
        reads.push_back("ACNCACACGATC"); // perfect across edge, pos 20
        reads.push_back("ACNCAACGATCG"); // 1 del across edge, pos 20
        reads.push_back("ACNCACCACGAT"); // 1 ins across edge, pos 20
        reads.push_back("ACTTGCTNCAGT"); // pos 1 insertion

        vargas::Aligner a(50, 12, 2, 2, 3, 1);
        auto res = a.align(reads, g.cbegin(), g.cend());
        REQUIRE(res.size() == 7);

        CHECK(res.max_score[0] == 22);
        CHECK(res.max_pos[0] == 12);
        CHECK(res.max_score[1] == 22);
        CHECK(res.max_pos[1] == 12);
        CHECK(res.max_score[2] == 15);
        CHECK(res.max_pos[2] == 16);
    }

}

TEST_CASE ("End to End alignment") {
    // Example from bowtie 2 manual
    vargas::Graph g;

    SUBCASE("BWT2 Local example") {
        /**
         * Read:      ACGGTTGCGTTAA-TCCGCCACG
         *                ||||||||| ||||||
         * Reference: TAACTTGCGTTAAATCCGCCTGG
         */
        const std::string read("ACGGTTGCGTTAATCCGCCACG"), ref("TAACTTGCGTTAAATCCGCCTGG");
        {
            vargas::Graph::Node n;
            n.set_as_ref();
            n.set_seq(ref);
            n.set_endpos(22); // 23 length -1 (for 0 indexed)
            g.add_node(n);
        }
        vargas::Aligner a(100, 22, 2, 6, 5, 3);
        auto res = a.align({read}, g.begin(), g.end());
        REQUIRE(res.size() == 1);
        CHECK(res.max_score[0] == 22);
        CHECK(res.max_pos[0] == 20);
    }

    SUBCASE("BWT2 ETE example") {
        /**
         * Read:      GACTGGGCGATCTCGACTTCG
         *            |||||  |||||||||| |||
         * Reference: GACTG--CGATCTCGACATCG
         */
        const std::string read("GACTGGGCGATCTCGACTTCG"), ref("GACTGCGATCTCGACATCG");
        {
            vargas::Graph::Node n;
            n.set_as_ref();
            n.set_seq(ref);
            n.set_endpos(18); // 19 length -1 (for 0 indexed) - 1 (pos of last base, not one after)
            g.add_node(n);
        }

        {
            vargas::AlignerETE a(100, 21, 0, 6, 5, 3);
            auto res = a.align({read}, g.begin(), g.end());
            REQUIRE(res.size() == 1);
            CHECK(res.max_pos[0] == 19);
            CHECK(res.max_score[0] == -17); // Best score -17 with bias 255
        }

        {
            vargas::WordAlignerETE a(100, 21, 0, 6, 5, 3);
            auto res = a.align({read}, g.begin(), g.end());
            REQUIRE(res.size() == 1);
            CHECK(res.max_pos[0] == 19);
            CHECK(res.max_score[0] == -17); // Best score -17 with bias 255
        }
    }

    SUBCASE("Bound check") {
        CHECK_THROWS(vargas::AlignerETE(100, 100, 3, 2, 2, 2));
    }
}

TEST_CASE ("Target score") {
    vargas::Graph g;
    vargas::Graph::Node n;
    n.set_seq("AAAACCCCCCCCCCCCAAA"); // Length 19
    n.set_endpos(18);
    g.add_node(n);

    const std::vector<std::string> reads = {"AAAA"};
    const std::vector<size_t> targets = {19};
    vargas::Aligner aligner(20, 4);
    auto res = aligner.align(reads, targets, g.begin(), g.end());
    REQUIRE(res.size() == 1);
    CHECK(res.max_score[0] == 8);
    CHECK(res.sub_score[0] == 6);
    CHECK(res.max_pos[0] == 4);
    CHECK(res.sub_pos[0] == 19);
    CHECK(res.correct[0] == 2);
    CHECK(res.target_score[0] == 6);
}

TEST_SUITE_END();

TEST_CASE ("prof") {
    simdpp::uint8x16 v;
    for (size_t i = 0; i < 16; ++i) vargas::insert<uint8_t>(rand() % 16, i, v);
    uint8_t prev = 0;
    auto start = std::chrono::steady_clock::now();
    for (size_t i = 0; i < 5000000000; ++i) {
        prev = vargas::extract<uint8_t>(prev, v);
    }
    auto stop = rg::chrono_duration(start);
    std::cerr << prev << stop << "s\n";
}
#endif //VARGAS_ALIGNMENT_H
