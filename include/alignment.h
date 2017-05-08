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
#include "simd.h"
#include "graph.h"
#include "doctest.h"
#include "simd.h"

#include <array>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <string>


#if VA_DEBUG_SW
#include <iomanip>
#endif

#define DEFAULT_TOL_FACTOR 4 // If the pos is +- read_len/tol, count as correct alignment

namespace vargas {


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
      virtual void set_scores(unsigned, unsigned, unsigned, unsigned) = 0;

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
      virtual void set_correctness_tolerance(const unsigned) = 0;

      /**
       * @return Current correctness tolerance
       */
      virtual unsigned tolerance() const = 0;

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
      virtual void align_into(const std::vector<std::string> &, const std::vector<unsigned> &,
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
      virtual Results align(const std::vector<std::string> &read_group, const std::vector<unsigned> &targets,
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
          std::vector<unsigned> targets(read_group.size());
          std::fill(targets.begin(), targets.end(), 0);
          return align(read_group, targets, begin, end);
      }

      static constexpr unsigned default_tolerance() { return DEFAULT_TOL_FACTOR; }

    protected:
      ScoreProfile _prof;

      /**
       * @brief
       * Ending vectors from a previous node
       */
      template<typename simd_t>
      struct _seed {
          _seed() = delete;
          _seed(const unsigned _read_len) : S_col(_read_len + 1), I_col(_read_len + 1) {}
          SIMDVector<simd_t> S_col; /**< Last column of score matrix.*/
          SIMDVector<simd_t> I_col;
      };

      template<typename simd_t>
      struct _target {
          unsigned char idx[simd_t::length + 1];
          int score[simd_t::length + 1];
          unsigned pos[simd_t::length + 1];
      };

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
   * Vargas::Aligner a(4);
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
   * @tparam native_t Native data type of score matrix element. One of uint8_t, uint16_t, uint32_t
   * @tparam END_TO_END If true, perform end to end alignment
   */
  template<typename simd_t, bool END_TO_END>
  class AlignerT: public AlignerBase {
    public:

      using native_t = typename simd_t::native_t;

      AlignerT(unsigned read_len, const ScoreProfile &prof) :
      _alignment_group(read_len),
      _S(read_len + 1), _Dc(read_len + 1), _Ic(read_len + 1),
      _read_len(read_len) {
          set_scores(prof); // May throw
          set_correctness_tolerance(read_len / DEFAULT_TOL_FACTOR);
      }

      /**
       * @brief
       * Set scoring parameters.
       * @param read_len max read length
       * @param match match score
       * @param mismatch mismatch penalty
       * @param open gap open penalty
       * @param extend gap extend penalty
       */
      AlignerT(unsigned read_len, unsigned match = 2, unsigned mismatch = 2, unsigned open = 3, unsigned extend = 1) :
      AlignerT(read_len, ScoreProfile(match, mismatch, open, extend)) {}


      template<typename S, bool E>
      AlignerT(const AlignerT<S, E> &a) = delete;

      template<typename S, bool E>
      AlignerT(AlignerT<S, E> &&a) = delete;

      template<typename S, bool E>
      AlignerT &operator=(const AlignerT<S, E> &) = delete;

      template<typename S, bool E>
      AlignerT &operator=(AlignerT<S, E> &&) = delete;

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

          AlignmentGroup(unsigned read_len) : _packaged_reads(read_len), _read_len(read_len) {}

          __RG_STRONG_INLINE__
          void load_reads(const std::vector<std::string> &reads, const unsigned begin, const unsigned end) {
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
           * Return the i'th base of every read in a simd vector.
           * @param i base index.
           */
          const simd_t &at(const unsigned i) const {
              return _packaged_reads.at(i);
          }

          /**
           * @brief
           * Pointer to raw packaged read data.
           * @return VEC_TYPE pointer
           */
          const simd_t *data() const {
              return _packaged_reads.data();
          }

          /**
           * @brief
           * Non const version of at(i).
           * @param i base index
           */
          simd_t &operator[](const int i) {
              return _packaged_reads[i];
          }

          const simd_t &operator[](const int i) const {
              return _packaged_reads[i];
          }

          /**
           * @brief
           * Returns optimal number of reads in a batch based on SIMD architecture.
           * @return batch size.
           */
          static constexpr unsigned group_size() { return simd_t::length; }

          /**
           * @return iterator to the beginning of the packaged reads.
           */
          typename SIMDVector<simd_t>::const_iterator begin() const {
              return _packaged_reads.cbegin();
          }

          /**
           * @return iterator to the end of the packaged reads.
           */
          typename SIMDVector<simd_t>::const_iterator end() const {
              return _packaged_reads.cend();
          }

        private:


          /**
           * _packaged_reads[i] contains all i'th bases.
           * The length of _packaged_reads is the length of the read,
           * where as the length of _packaged_reads[i] is the number
           * of reads.
           */
          SIMDVector<simd_t> _packaged_reads;
          const unsigned _read_len;

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
              for (unsigned r = 0; r < _reads.size(); ++r) {
                  assert(_reads[r].size() == _read_len);
                  // Put each base in the appropriate vector element
                  for (unsigned p = 0; p < _read_len; ++p) {
                      _packaged_reads[p][r] = _reads[r][p];
                  }
              }

              // Pad underful batches
              for (unsigned r = _reads.size(); r < group_size(); ++r) {
                  for (unsigned p = 0; p < _read_len; ++p) {
                      _packaged_reads[p][r] = rg::Base::N;
                  }
              }
          }

      };

      void set_scores(unsigned match, unsigned mismatch, unsigned open, unsigned extend) override {
          _prof = ScoreProfile(match, mismatch, open, extend);
          set_scores(_prof);
      }

      void set_scores(const ScoreProfile &prof) override {
          if (!END_TO_END && (prof.match < 0 || prof.read_gopen < 0 ||
          prof.read_gext < 0 || prof.ref_gopen < 0 || prof.ref_gext < 0)) {
              throw std::invalid_argument("Expected match bonus and penalties > 0 in local alignment mode.");
          }
          _prof = prof;
          _prof.end_to_end = END_TO_END;
          _match_vec = prof.match;
          _mismatch_vec = -prof.mismatch;
          _gap_open_extend_vec_rd = prof.read_gopen + prof.read_gext;
          _gap_extend_vec_rd = prof.read_gext;
          _gap_open_extend_vec_ref = prof.ref_gopen + prof.ref_gext;
          _gap_extend_vec_ref = prof.ref_gext;
          _ambig_vec = -prof.ambig;
          _bias = _get_bias(_read_len, prof.match, prof.mismatch, prof.read_gopen, prof.read_gext);
          _Dc[0] = _bias;
          set_correctness_tolerance(prof.tol);
      }

      void set_correctness_tolerance(const unsigned tol) override {
          _prof.tol = tol;
      }

      unsigned tolerance() const override {
          return _prof.tol;
      }

      /**
       * @return maximum number of reads that can be aligned at once.
       */
      static constexpr unsigned read_capacity() { return simd_t::length; }

      void align_into(const std::vector<std::string> &read_group, const std::vector<pos_t> &targets,
                      Graph::const_iterator begin, Graph::const_iterator end, Results &aligns) override {



          assert(targets.size() == read_group.size());

          const unsigned num_groups = 1 + ((read_group.size() - 1) / read_capacity());
          // Possible oversize if there is a partial group
          aligns.resize(num_groups * read_capacity());

          // Keep the scores at the positions, overwrites position. [0] is current position, 1-:ead_capacity + 1 is pos
          std::unordered_map<unsigned, _seed<simd_t>> seed_map; // Maps node ID to the ending matrix cols of the node
          _seed <simd_t> seed(_read_len);
          const std::vector<unsigned> zvec = {0};

          for (unsigned group = 0; group < num_groups; ++group) {
              seed_map.clear();

              // Subset of read set
              const unsigned beg_offset = group * read_capacity();
              const unsigned end_offset = std::min<unsigned>((group + 1) * read_capacity(), read_group.size());
              const unsigned len = end_offset - beg_offset;
              assert(len <= read_capacity());

              _alignment_group.load_reads(read_group, beg_offset, end_offset);
              _max_score = std::numeric_limits<native_t>::min();
              _max_pos = aligns.max_pos.data() + beg_offset;

              _sub_score = std::numeric_limits<native_t>::min();
              _sub_pos = aligns.sub_pos.data() + beg_offset;
              _max_count = aligns.max_count.data() + beg_offset;
              _sub_count = aligns.sub_count.data() + beg_offset;

              // Create subrange of targets, and sort by position.
              for (unsigned j = 0; j < len; ++j) {
                  _target_subrange.idx[j] = j;
                  _target_subrange.pos[j] = targets[beg_offset + j];
                  _target_subrange.score[j] = std::numeric_limits<int>::min();
              }
              // pad with extra
              for (unsigned j = len; j < read_capacity() + 1; ++j) {
                  _target_subrange.pos[j] = std::numeric_limits<int>::max();
              }

              // Sort ascending position
              for (unsigned c = 0; c < len - 1; ++c) {
                  for (unsigned d = 0; d < len - c - 1; ++d) {
                      if (_target_subrange.pos[d] > _target_subrange.pos[d + 1]) {
                          std::swap(_target_subrange.pos[d], _target_subrange.pos[d + 1]);
                          std::swap(_target_subrange.idx[d], _target_subrange.idx[d + 1]);
                      }
                  }
              }

              for (auto gi = begin; gi != end; ++gi) {
                  _get_seed(gi.incoming(), seed_map, seed);
                  if (gi->is_pinched()) seed_map.clear();
                  _fill_node(*gi, _alignment_group, seed, seed_map.emplace(gi->id(), _read_len).first->second);
              }



              // Copy scores
              for (unsigned char i = 0; i < len; ++i) {
                  aligns.max_score[beg_offset + i] = _max_score[i] - _bias;
                  aligns.sub_score[beg_offset + i] = _sub_score[i] - _bias;
                  aligns.target_score[beg_offset + _target_subrange.idx[i]] = _target_subrange.score[i] - _bias;
              }

          }
          // Crop off potential buffer
          aligns.resize(read_group.size());
          aligns.profile = _prof;
          aligns.finalize(targets);

      }


    private:

      void _seed_matrix(_seed <simd_t> &seed) const {
          if (END_TO_END) {
              for (unsigned i = 0; i < _read_len; ++i) {
                  int v = _bias - _prof.ref_gopen - ((i + 1) * _prof.ref_gext);
                  seed.S_col[i + 1] =
                  v < std::numeric_limits<native_t>::min() ? std::numeric_limits<native_t>::min() : v;
              }
          } else {
              std::fill(seed.S_col.begin(), seed.S_col.end(), _bias);
          }
          seed.I_col = seed.S_col;
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
      void _get_seed(const std::vector<unsigned> &prev_ids,
                     std::unordered_map<unsigned, _seed<simd_t>> &seed_map, _seed<simd_t> &seed) const {
          if (prev_ids.size() == 0) {
              _seed_matrix(seed);
          }
          else {
              for (unsigned i = 1; i < _read_len + 1; ++i) {
                  const auto &s = seed_map.at(prev_ids[0]);
                  seed.S_col[i] = s.S_col[i];
                  seed.I_col[i] = s.I_col[i];

                  for (unsigned p = 1; p < prev_ids.size(); ++p) {
                      const auto &s = seed_map.at(prev_ids[p]);
                      seed.S_col[i] = max(seed.S_col[i], s.S_col[i]);
                      seed.I_col[i] = max(seed.I_col[i], s.I_col[i]);
                  }
              }
          }

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
      void _fill_node(const Graph::Node &n, const AlignmentGroup &read_group,
                      const _seed <simd_t> &s, _seed <simd_t> &nxt) {
          // Empty nodes represents deletions
          if (n.seq().size() == 0) {
              nxt = s;
              return;
          }

          unsigned curr_pos = n.end_pos() - n.seq().size() + 2;

          int csp = 0;
          while (_target_subrange.pos[csp] < curr_pos) ++csp;

          _S = s.S_col;
          _Ic = s.I_col;
          for (const rg::Base ref_base : n) {

              _Sd = _bias;
              for (unsigned r = 0; r < _read_len; ++r) {
                  _fill_cell(read_group[r], ref_base, r + 1, curr_pos);
              }
              if (END_TO_END) _fill_cell_finish(_read_len, curr_pos);

              while (_target_subrange.pos[csp] == curr_pos) {
                  for (unsigned q = END_TO_END ? _read_len : 1; q < _read_len + 1; ++q) {
                      _target_subrange.score[csp] = std::max<int>(_target_subrange.score[csp],
                                                                  _S[q][_target_subrange.idx[csp]]);
                  }
                  ++csp;
              }
              ++curr_pos;
          }

          nxt.S_col = _S;
          nxt.I_col = _Ic;

      }

      /**
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row _curr_posent row in matrix
       * @param col _curr_posent column in matrix
       */
      __RG_STRONG_INLINE__
      void _fill_cell(const simd_t &read, const rg::Base &ref, const unsigned &row, const pos_t &curr_pos) {
          _Dc[row] = max(_Dc[row - 1] - _gap_extend_vec_ref, _S[row - 1] - _gap_open_extend_vec_ref);
          _Ic[row] = max(_Ic[row] - _gap_extend_vec_rd, _S[row] - _gap_open_extend_vec_rd);
          simd_t sr;
          if (ref != rg::Base::N) {
              sr = blend(read == ref, _match_vec, _mismatch_vec);
              sr = _Sd + blend(read == rg::Base::N, _ambig_vec, sr);
          } else {
              sr = _Sd + _ambig_vec;
          }

          _Sd = _S[row]; // S(i-1, j-1) for the next cell to be filled in
          _S[row] = max(_Ic[row], max(_Dc[row], sr));

          if (!END_TO_END) _fill_cell_finish(row, curr_pos);
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
      void _fill_cell_finish(const unsigned &row, const pos_t &curr_pos) {
          simd_t _tmp0;
          _tmp0 = _S[row] == _max_score;
          if (_tmp0) {
              // Check for equal max score.
              for (unsigned i = 0; i < read_capacity(); ++i) {
                  if (_tmp0[i]) {
                      if (curr_pos > _max_pos[i] + _read_len) ++(_max_count[i]);
                      _max_pos[i] = curr_pos;
                  }
              }
          }

          _tmp0 = _S[row] > _max_score;
          if (_tmp0) {
              // Check for new or equal high scores
              _max_score = max(_S[row], _max_score);
              for (unsigned i = 0; i < read_capacity(); ++i) {
                  if (_tmp0[i]) {
                      // Demote old max to submax
                      if (curr_pos > _max_pos[i] + _read_len) {
                          _sub_score[i] = _max_score[i];
                          _sub_pos[i] = _max_pos[i];
                          _sub_count[i] = _max_count[i];
                      }
                      _max_count[i] = 1;
                      _max_pos[i] = curr_pos;
                  }
              }
          }


          _tmp0 = _S[row] == _sub_score;
          if (_tmp0) {
              // Repeat sub score
              for (unsigned i = 0; i < read_capacity(); ++i) {
                  if (_tmp0[i] && curr_pos > _max_pos[i] + _read_len) {
                      _sub_count[i] += curr_pos > (_sub_pos[i] + _read_len);
                      _sub_pos[i] = curr_pos;
                  }
              }
          }

          // Greater than old sub max and less than max score (prevent repeats of max triggering)
          _tmp0 = (_S[row] > _sub_score) & (_S[row] < _max_score);
          if (_tmp0) {
              // new second best score
              for (unsigned i = 0; i < read_capacity(); ++i) {
                  if (_tmp0[i] && curr_pos > _max_pos[i] + _read_len) {
                      _sub_score[i] = _S[row][i];
                      _sub_count[i] = 1;
                      _sub_pos[i] = curr_pos;
                  }
              }
          }


      }

      static native_t _get_bias(const unsigned read_len, const unsigned match, const unsigned mismatch,
                                const unsigned gopen, const unsigned gext) {
          static bool has_warned = false;
          if (read_len * match > std::numeric_limits<native_t>::max() - std::numeric_limits<native_t>::min()) {
              throw std::domain_error("Insufficient bit-width for given match score and read length.");
          }
          if (!END_TO_END) return std::numeric_limits<native_t>::min();

          // End to end alignment
          unsigned int b = std::numeric_limits<native_t>::max() - (read_len * match);

          //TODO Could be relaxed - all indels or all mismatch is unreasonable
          if (!has_warned && (gopen + (gext * (read_len - 1)) > b || read_len * mismatch > b)) {
              std::cerr << "WARN: Possibility of score saturation with parameters in end-to-end mode:\n\t"
                        << "Cell Width: "
                        << (int) std::numeric_limits<native_t>::max() - (int) std::numeric_limits<native_t>::min()
                        << ", Bias: " << b << "\n";
              has_warned = true;
          }
          return b;
      }


      /*********************************** Variables ***********************************/

      AlignmentGroup _alignment_group;
      SIMDVector<simd_t> _S, _Dc, _Ic;
      _target<simd_t> _target_subrange; // Pad with one max so serve as buffer

      simd_t
      _match_vec, _mismatch_vec, _ambig_vec,
      _gap_open_extend_vec_rd, _gap_extend_vec_rd,
      _gap_open_extend_vec_ref, _gap_extend_vec_ref,
      _Sd, _max_score, _sub_score;

      pos_t *_max_pos, *_sub_pos;
      unsigned  *_max_count, *_sub_count;

      native_t _bias;
      const unsigned int _read_len;

  };

  using Aligner = AlignerT<int8_fast, false>;
  using WordAligner = AlignerT<int16_fast, false>;
  using AlignerETE = AlignerT<int8_fast, true>;
  using WordAlignerETE = AlignerT<int16_fast, true>;

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
        const std::vector<unsigned> origins = {8, 8, 5, 5, 7, 6, 10, 6};

        vargas::Results aligns;
        {
            vargas::Aligner a(7);
            aligns = a.align(reads, origins, g.begin(), g.end());
        }
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
        CHECK((int) aligns.correct[7] == 0); // lower cost to end after AAAG
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
        const std::vector<unsigned> origins = {8, 8, 5, 5, 7, 6, 10, 4, 10, 10};

        // hisat like params
        vargas::Aligner a(10, 2, 6, 5, 3);
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
        vargas::Aligner a(10, prof);
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
        const std::vector<unsigned> origins = {8, 8, 5, 5, 7, 6, 10, 6};

        vargas::WordAligner a(7);
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
        CHECK((int) aligns.correct[7] == 0);
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
        const std::vector<unsigned> origins = {8, 8, 5, 5, 7, 6, 10, 4, 10, 10};

        // hisat like params
        vargas::WordAligner a(10, 2, 6, 5, 3);
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
        n.set_endpos(67);
        n.set_as_ref();
        n.set_seq("ACGATCGTACGCNAGCTAGCCACAGTGCCCCCCTATATACGAN");
        g.add_node(n);
    }
    g.add_edge(0, 1);

    SUBCASE ("Indel") {
        std::vector<std::string> reads;
        reads.push_back("ACTGCTNCAGTC"); // perfect alignment, pos 1
        reads.push_back("ACTGCTACAGTC"); // perfect alignment, pos 1, diff N
        reads.push_back("CCACAGCCCCCC"); // 2 del
        reads.push_back("ACNCACACGATC"); // perfect across edge
        reads.push_back("ACNCAACGATCG"); // 1 del across edge
        reads.push_back("ACNCACCACGAT"); // 1 ins across edge
        reads.push_back("ACTTGCTNCAGT"); // 1 ins
        reads.push_back("ACNCACCGATCG");
        reads.push_back("NACNCAACGATC");
        reads.push_back("AGCCTTACAGTG"); // 2 ins

        SUBCASE("Same read/ref") {
            vargas::Aligner a(12, 2, 6, 3, 1);
            auto res = a.align(reads, g.begin(), g.end());
            REQUIRE(res.size() == 10);

            CHECK(res.max_score[0] == 22);
            CHECK(res.max_pos[0] == 12);
            CHECK(res.max_score[1] == 22);
            CHECK(res.max_pos[1] == 12);
            CHECK(res.max_score[2] == 19);
            CHECK(res.max_pos[2] == 58);
            CHECK(res.max_score[3] == 22);
            CHECK(res.max_pos[3] == 31);
            CHECK(res.max_score[4] == 18);
            CHECK(res.max_pos[4] == 32);
            CHECK(res.max_score[5] == 16);
            CHECK(res.max_pos[5] == 30);
            CHECK(res.max_score[6] == 16);
            CHECK(res.max_pos[6] == 11);
            CHECK(res.max_score[7] == 18);
            CHECK(res.max_pos[7] == 32);
            CHECK(res.max_score[8] == 16);
            CHECK(res.max_pos[8] == 31);
            CHECK(res.max_score[9] == 15);
            CHECK(res.max_pos[9] == 52);

        }

        SUBCASE("Diff read/ref") {
            vargas::ScoreProfile prof(2, 6, 4, 1, 2, 1);
            vargas::Aligner a(12, prof);
            auto res = a.align(reads, g.begin(), g.end());
            REQUIRE(res.size() == 10);

            CHECK(res.max_score[0] == 22);
            CHECK(res.max_pos[0] == 12);
            CHECK(res.max_score[1] == 22);
            CHECK(res.max_pos[1] == 12);
            CHECK(res.max_score[2] == 18);
            CHECK(res.max_pos[2] == 58);
            CHECK(res.max_score[3] == 22);
            CHECK(res.max_pos[3] == 31);
            CHECK(res.max_score[4] == 17);
            CHECK(res.max_pos[4] == 32);
            CHECK(res.max_score[5] == 17);
            CHECK(res.max_pos[5] == 30);
            CHECK(res.max_score[6] == 17);
            CHECK(res.max_pos[6] == 11);
            CHECK(res.max_score[7] == 17);
            CHECK(res.max_pos[7] == 32);
            CHECK(res.max_score[8] == 15);
            CHECK(res.max_pos[8] == 31);
            CHECK(res.max_score[9] == 16);
            CHECK(res.max_pos[9] == 52);
        }
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
        vargas::Aligner a(22, 2, 6, 5, 3);
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
            vargas::AlignerETE a(21, 0, 6, 5, 3);
            auto res = a.align({read}, g.begin(), g.end());
            REQUIRE(res.size() == 1);
            CHECK(res.max_pos[0] == 19);
            CHECK(res.max_score[0] == -17); // Best score -17 with bias 255
        }

        {
            vargas::WordAlignerETE a(21, 0, 6, 5, 3);
            auto res = a.align({read}, g.begin(), g.end());
            REQUIRE(res.size() == 1);
            CHECK(res.max_pos[0] == 19);
            CHECK(res.max_score[0] == -17); // Best score -17 with bias 255
        }
    }

    SUBCASE("Bound check") {
        CHECK_THROWS(vargas::AlignerETE(100, 3, 2, 2, 2));
    }
}

TEST_CASE ("Target score") {
    vargas::Graph g;
    vargas::Graph::Node n;
    n.set_seq("AAAACCCCCCCCCCCCAAA"); // Length 19
    n.set_endpos(18);
    g.add_node(n);

    const std::vector<std::string> reads = {"AAAA"};
    const std::vector<unsigned> targets = {19};
    vargas::Aligner aligner(4);
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

#endif //VARGAS_ALIGNMENT_H
