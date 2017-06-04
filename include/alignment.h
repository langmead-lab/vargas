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
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
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
#include <stdexcept>

/**
 * Any debug level prints SW grid for VARGAS_ALIGN_DEBUG_N element of SIMD vector.
 * > 1 prints query profile.
 */
#define VARGAS_ALIGN_DEBUG 0
#define VARGAS_ALIGN_DEBUG_N 0 // Which SIMD element

namespace vargas {

  /**
   * @brief
   * Common base class to all template versions of Aligners
   */
  class AlignerBase {
    public:
      virtual ~AlignerBase() = 0;

      /**
       * @brief
       * Set the score profile from a ScoreProfile - gap penalties between ref/read may vary.
       * @param prof
       */
      virtual void set_scores(const ScoreProfile &prof) = 0;

      /**
       * @brief
       * Align a batch of reads to a graph range, return a vector of alignments
       * corresponding to the reads.
       * @param read_group vector of reads to align to
       * @param quals Quality values
       * @param begin iterator to beginning of graph
       * @param end iterator to end of graph
       * @param aligns Results packet to populate
       * @param fwdonly Only align to forward strand
       */
      virtual void align_into(const std::vector<std::string> &,
                              const std::vector<std::vector<char>> &,
                              Graph::const_iterator, Graph::const_iterator, Results &, bool) = 0;

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
                            Graph::const_iterator begin, Graph::const_iterator end, bool fwdonly=true) {
          Results aligns;
          align_into(read_group, {}, begin, end, aligns, fwdonly);
          return aligns;
      }

    protected:
      ScoreProfile _prof;

      /**
       * @brief
       * Ending vectors from a previous node
       */
      template<typename st>
      struct _seed {
          _seed() = delete;
          _seed(const unsigned _read_len) : S_col(_read_len + 1), I_col(_read_len + 1) {}
          SIMDVector<st> S_col; /**< Last column of score matrix.*/
          SIMDVector<st> I_col;
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
   * @tparam simd_t data type of score matrix element. One of SIMD<uint8_t>, SIMD<uint16_t>
   * @tparam END_TO_END If true, perform end to end alignment
   * @tparam MSONLY Only collect max score- no positions or subscores
   */
  template<typename simd_t, bool END_TO_END, bool MSONLY=false>
  class AlignerT: public AlignerBase {
    public:

      using native_t = typename simd_t::native_t;
      /**
       * @brief
       * Query profile type
       * @details
       * Index levels: \n
       * 0 : indicates base of read
       * 1 : Reference base, one of A,C,G,T,N
       * 2 : mismatch penalties for each read
       */
      using qp_t = std::vector<std::array<simd_t, 5>, aligned_allocator<std::array<simd_t, 5>, simd_t::size>>;

      AlignerT(unsigned read_len, const ScoreProfile &prof) :
      _alignment_group(read_len),
      _S(read_len + 1), _Dc(read_len + 1), _Ic(read_len + 1),
      _read_len(read_len) {
          set_scores(prof); // May throw
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
      AlignerT(unsigned read_len,unsigned match = 2, unsigned mismatch = 2, unsigned open = 3, unsigned extend = 1) :
      AlignerT(read_len, ScoreProfile(match, mismatch, open, extend)) {}


      template<typename S, bool E, bool M> AlignerT(const AlignerT<S, E, M> &a) = delete;
      template<typename S, bool E, bool M> AlignerT(AlignerT<S, E, M> &&a) = delete;
      template<typename S, bool E, bool M> AlignerT &operator=(const AlignerT<S, E, M> &) = delete;
      template<typename S, bool E, bool M> AlignerT &operator=(AlignerT<S, E, M> &&) = delete;

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

          AlignmentGroup(unsigned read_len) : _query_prof(read_len), _rd_ln(read_len) {}
          AlignmentGroup() = delete;

          void load_reads(const std::vector<std::string> &reads,
                          const std::vector<std::vector<char>> &quals,
                          const ScoreProfile &prof,
                          const unsigned begin, const unsigned end, bool revcomp = false) {
              load_reads(std::vector<std::string>(reads.begin() + begin, reads.begin() + end), quals, prof, revcomp);
          }

          /**
           * @param batch load the given vector of reads.
           */
          void load_reads(const std::vector<std::string> &batch,
                          const std::vector<std::vector<char>> &quals,
                          const ScoreProfile &prof,
                          bool revcomp) {
              std::vector<std::vector<rg::Base>> _reads;
              for (auto &b : batch) _reads.push_back(rg::seq_to_num(b));
              load_reads(_reads, quals, prof, revcomp);
          }

          /**
           * @param batch load the given vector of reads.
           */
          void load_reads(const std::vector<std::vector<rg::Base>> &batch,
                          const std::vector<std::vector<char>> &quals,
                          const ScoreProfile &prof,
                          bool revcomp) {
              _package_reads(batch, quals, prof, revcomp);
          }

          typename qp_t::value_type &at(const unsigned i) const {
              return _query_prof.at(i);
          }

          typename qp_t::value_type &operator[](const int i) {
              return _query_prof[i];
          }

          const typename qp_t::value_type &operator[](const int i) const {
              return _query_prof[i];
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
          typename qp_t::const_iterator begin() const {
              return _query_prof.cbegin();
          }

          /**
           * @return iterator to the end of the packaged reads.
           */
          typename qp_t::const_iterator end() const {
              return _query_prof.cend();
          }

          const qp_t &query_profile() const {
              return _query_prof;
          }

        private:
          qp_t _query_prof;
          const unsigned _rd_ln;

          /**
           * Interleaves reads so all same-index base positions are in one
           * vector. Empty spaces are padded with Base::N.
           * @param reads vector of reads to package
           * @param quals Quality values for reads, Phred
           * @param prof ScoreProfile
           * @param revcomp Use reverse complement
           */
          void _package_reads(const std::vector<std::vector<rg::Base>> &reads,
                              const std::vector<std::vector<char>> &quals,
                              const ScoreProfile &prof,
                              bool revcomp) {
              assert(reads.size() <= group_size());
              static constexpr std::array<rg::Base, 4> bases = {rg::Base::A, rg::Base::C, rg::Base::G, rg::Base::T};
              // Interleave reads
              // For each read (read[i] is in _packaged_reads[0..n][i]

              const int inc = revcomp ? -1 : 1;

              for (unsigned r = 0; r < reads.size(); ++r) {

                  // Prepend short reads with 0
                  int pos = _rd_ln - reads[r].size();
                  for (int i = 0; i < pos; ++i) {
                      for (auto b : bases) _query_prof[i][b][r] = 0;
                  }

                  const int start = revcomp ? reads[r].size() - 1 : 0;
                  const int end = revcomp ? -1 : reads[r].size();

                  // Store in index pos, run through bases from idx start --> end based on revcomp
                  for (int p = start; p != end; p += inc) {
                      const auto rdb = revcomp ? rg::complement_b(reads[r][p]) : reads[r][p];
                      _query_prof[pos][rg::Base::N][r] = -prof.ambig;
                      for (auto b : bases) {
                          if (rdb == rg::Base::N) _query_prof[pos][b][r] = -prof.ambig;
                          else if (rdb == b) _query_prof[pos][b][r] = prof.match;
                          else if (!quals.size() || !quals[r].size()) _query_prof[pos][b][r] = -prof.mismatch_max;
                          else _query_prof[pos][b][r] = -prof.penalty(quals[r][p]);
                      }
                      ++pos;
                  }
              }

              #if VARGAS_ALIGN_DEBUG > 1
              // Print query profile
              std::cerr << "Query Profile (" << VARGAS_ALIGN_DEBUG_N << ")" << std::endl;
              for (unsigned i = 0; i < _rd_ln - reads[VARGAS_ALIGN_DEBUG_N].size(); ++i) std::cerr << "\t-";
              for (auto b : reads[VARGAS_ALIGN_DEBUG_N]) std::cerr << "\t" << rg::num_to_base(b);
              std::cerr << std::endl << "A\t";
              for (unsigned i = 0; i < _rd_ln; ++i) {
                  std::cerr << (int)_query_prof[i][rg::Base::A][VARGAS_ALIGN_DEBUG_N] << '\t';
              }
              std::cerr << std::endl << "C\t";

              for (unsigned i = 0; i < _rd_ln; ++i) {
                  std::cerr << (int)_query_prof[i][rg::Base::C][VARGAS_ALIGN_DEBUG_N] << '\t';
              }
              std::cerr << std::endl << "G\t";

              for (unsigned i = 0; i < _rd_ln; ++i) {
                  std::cerr << (int)_query_prof[i][rg::Base::G][VARGAS_ALIGN_DEBUG_N] << '\t';
              }
              std::cerr << std::endl << "T\t";

              for (unsigned i = 0; i < _rd_ln; ++i) {
                  std::cerr << (int)_query_prof[i][rg::Base::T][VARGAS_ALIGN_DEBUG_N] << '\t';
              }
              std::cerr << std::endl << "N\t";

              for (unsigned i = 0; i < _rd_ln; ++i) {
                  std::cerr << (int)_query_prof[i][rg::Base::N][VARGAS_ALIGN_DEBUG_N] << '\t';
              }
              std::cerr << std::endl << std::endl;
              #endif
          }

      };


      void set_scores(const ScoreProfile &prof) override {
          _prof = prof;
          _prof.end_to_end = END_TO_END;
          _bias = _get_bias(_read_len, prof.match, prof.mismatch_max, prof.read_gopen, prof.read_gext);
          _Dc[0] = _S[0] = _bias;
          _gap_extend_vec_rd = prof.read_gext;
          _gap_extend_vec_ref = prof.ref_gext;
          _gap_open_extend_vec_rd = prof.read_gopen + prof.read_gext;
          _gap_open_extend_vec_ref = prof.ref_gopen + prof.ref_gext;
      }

      /**
       * @return maximum number of reads that can be aligned at once.
       */
      static constexpr unsigned read_capacity() { return simd_t::length; }

      void align_into(const std::vector<std::string> &read_group,
                      const std::vector<std::vector<char>> &quals,
                      Graph::const_iterator begin, Graph::const_iterator end,
                      Results &aligns, bool fwdonly=true) override {

          const unsigned num_groups = 1 + ((read_group.size() - 1) / read_capacity());
          // Possible oversize if there is a partial group
          aligns.resize(num_groups * read_capacity());

          // Keep the scores at the positions, overwrites position. [0] is current position, 1-:ead_capacity + 1 is pos
          std::unordered_map<unsigned, _seed<simd_t>> seed_map; // Maps node ID to the ending matrix cols of the node
          _seed <simd_t> seed(_read_len);

          if (fwdonly) {
              std::fill(aligns.max_strand.begin(), aligns.max_strand.end(), Strand::FWD);
              std::fill(aligns.sub_strand.begin(), aligns.sub_strand.end(), Strand::FWD);
          }

          for (unsigned group = 0; group < num_groups; ++group) {
              seed_map.clear();

              // Subset of read set
              const unsigned beg_offset = group * read_capacity();
              const unsigned end_offset = std::min<unsigned>((group + 1) * read_capacity(), read_group.size());
              const unsigned len = end_offset - beg_offset;
              assert(len <= read_capacity());

              _max_score = std::numeric_limits<native_t>::min();

              if (!MSONLY) {
                  _max_pos = aligns.max_pos.data() + beg_offset;
                  _sub_score = std::numeric_limits<native_t>::min();
                  _sub_pos = aligns.sub_pos.data() + beg_offset;
                  _max_count = aligns.max_count.data() + beg_offset;
                  _sub_count = aligns.sub_count.data() + beg_offset;
              }

              // Forward
              _alignment_group.load_reads(read_group, quals, _prof, beg_offset, end_offset, false);
              for (auto gi = begin; gi != end; ++gi) {
                  _get_seed(gi.incoming(), seed_map, seed);
                  if (gi->is_pinched()) seed_map.clear();
                  _fill_node(*gi, _alignment_group.query_profile(), seed, seed_map.emplace(gi->id(), _read_len).first->second);
              }

              // Reverse
              if (!fwdonly) {
                  seed_map.clear();
                  _alignment_group.load_reads(read_group, quals, _prof, beg_offset, end_offset, true);
                  simd_t fwdmax = _max_score;
                  simd_t fwdsub = _sub_score;

                  for (auto gi = begin; gi != end; ++gi) {
                      _get_seed(gi.incoming(), seed_map, seed);
                      if (gi->is_pinched()) seed_map.clear();
                      _fill_node(*gi, _alignment_group.query_profile(), seed, seed_map.emplace(gi->id(), _read_len).first->second);
                  }

                  // Assign strands
                  for(size_t i = 0; i < len; ++i) {
                     aligns.max_strand[beg_offset + i] = _max_score[i] > fwdmax[i] ? Strand::REV : Strand::FWD;
                     aligns.sub_strand[beg_offset + i] = _sub_score[i] > fwdsub[i] ? Strand::REV : Strand::FWD;
                  }
              }


              // Copy scores
              for (unsigned char i = 0; i < len; ++i) {
                  aligns.max_score[beg_offset + i] = _max_score[i] - _bias;
                  if (!MSONLY) {
                      aligns.sub_score[beg_offset + i] = _sub_score[i] - _bias;
                  }
              }

          }
          // Crop off potential buffer
          aligns.resize(read_group.size());
          aligns.profile = _prof;

      }

    private:

      /**
       * @brief
       * Seeds the matrix when there are no previous nodes. In end to end mode, the seed is penalized.
       * Otherwise, seeds are cleared to _bias.
       * @param seed
       */
      void _seed_matrix(_seed <simd_t> &seed) const {
          if (END_TO_END) {
              seed.S_col[0] = _bias;
              for (unsigned i = 1; i <= _read_len; ++i) {
                  int v = _bias - _prof.ref_gopen - (i * _prof.ref_gext);
                  seed.S_col[i] = v < std::numeric_limits<native_t>::min() ? std::numeric_limits<native_t>::min() : v;
              }
          }
          else {
              std::fill(seed.S_col.begin(), seed.S_col.end(), _bias);
          }
          seed.I_col = seed.S_col;
      }

      /**
       * @brief
       * Returns the best seed from all previous nodes.
       * Graph should be validated before alignment to ensure proper seed fetch
       * @param prev_ids All nodes preceding _curr_pos node
       * @param seed_map ID->seed map for all previous nodes
       * @param seed best seed to populate
       * @throws std::domain_error if a node listed as a previous node but it has not been encountered yet.
       * i.e. not topographically sorted.
       */
      __RG_STRONG_INLINE__
      void _get_seed(const std::vector<unsigned> &prev_ids, std::unordered_map<unsigned, _seed<simd_t>> &seed_map,
                     _seed<simd_t> &seed) const {
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
      void _fill_node(const Graph::Node &n, const qp_t &read_group,
                      const _seed <simd_t> &s, _seed <simd_t> &nxt) {
          // Empty nodes represents deletions
          if (n.seq().size() == 0) {
              nxt = s;
              return;
          }

          #if VARGAS_ALIGN_DEBUG
          std::vector<std::vector<char>> grid;
          grid.resize(_read_len);
          for (auto &&r : grid) r.resize(n.seq().size());
          unsigned deb_col = 0;
          #endif

          unsigned curr_pos = n.end_pos() - n.seq().size() + 2;

          _S = s.S_col;
          _Ic = s.I_col;
          for (const rg::Base ref_base : n) {
              _Sd = _bias;

              for (unsigned r = 0; r < _read_len; ++r) {
                  _fill_cell(read_group[r], ref_base, r + 1, curr_pos);
                  #if VARGAS_ALIGN_DEBUG
                  grid[r][deb_col] = _S[r+1][VARGAS_ALIGN_DEBUG_N];
                  #endif
              }
              if (END_TO_END) _fill_cell_finish(_read_len, curr_pos);
              ++curr_pos;

              #if VARGAS_ALIGN_DEBUG
              ++deb_col;
              #endif
          }

          #if VARGAS_ALIGN_DEBUG
          std::cerr << std::endl << "S";
          for (auto b : n) std::cerr << '\t' << rg::num_to_base(b);
          std::cerr << std::endl;
          for (unsigned i = 0; i < grid.size(); ++i) {
              const auto &row = grid[i];
              auto tmp = s; // s is const, make copy
              std::cerr << (int) tmp.S_col.at(i)[VARGAS_ALIGN_DEBUG_N];
              for (auto &&b : row) {
                  std::cerr << '\t' << (int) b;
              }
              std::cerr << '\n';
          }
          #endif

          nxt.S_col = _S;
          nxt.I_col = _Ic;
      }

      /**
       * @param read_base ReadBatch vector
       * @param ref reference sequence base
       * @param row _curr_pos row in matrix
       * @param col _curr_pos column in matrix
       */
      __RG_STRONG_INLINE__
      void _fill_cell(const typename qp_t::value_type &prof, const rg::Base &ref,
                      const unsigned &row, const pos_t &curr_pos) {
          _Dc[row] = max(_Dc[row - 1] - _gap_extend_vec_ref, _S[row - 1] - _gap_open_extend_vec_ref);
          _Ic[row] = max(_Ic[row] - _gap_extend_vec_rd, _S[row] - _gap_open_extend_vec_rd);
          simd_t sr = _Sd + prof[ref];
          _Sd = _S[row]; // S(i-1, j-1) for the next cell to be filled in
          _S[row] = max(_Ic[row], max(_Dc[row], sr));
          if (!END_TO_END) _fill_cell_finish(row, curr_pos);
      }

      /**
       * @brief
       * Takes the max of D,I, and M vectors and stores the _curr_pos best score/position
       * Currently does not support non-default template args
       * @param row _curr_pos row
       * @param col _curr_pos column
       * @param node_origin Current position, used to get absolute alignment position
       */
      __RG_STRONG_INLINE__ __RG_UNROLL__
      void _fill_cell_finish(const unsigned &row, const pos_t &curr_pos) {
          if (MSONLY) {
              _max_score = max(_S[row], _max_score);
          }
          else {
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
              std::cerr << "[warn] Possibility of score saturation with parameters in end-to-end mode:\n"
                        << "\tCell Width: "
                        << (int) std::numeric_limits<native_t>::max() - (int) std::numeric_limits<native_t>::min()
                        << ", Bias: " << b << ", Limits: gaplen=" << (b - gopen)/gext << " OR mismatches=" << b/mismatch << "\n";
              has_warned = true;
          }
          return b;
      }

      /*********************************** Variables ***********************************/

      AlignmentGroup _alignment_group;
      SIMDVector<simd_t> _S, _Dc, _Ic;

      simd_t _Sd, _max_score, _sub_score,
      _gap_extend_vec_ref, _gap_open_extend_vec_ref, _gap_extend_vec_rd, _gap_open_extend_vec_rd;

      pos_t *_max_pos, *_sub_pos;
      unsigned  *_max_count, *_sub_count;

      native_t _bias;
      const unsigned int _read_len;

  };

  using Aligner = AlignerT<int8_fast, false, false>;
  using WordAligner = AlignerT<int16_fast, false, false>;
  using AlignerETE = AlignerT<int8_fast, true, false>;
  using WordAlignerETE = AlignerT<int16_fast, true, false>;

  using MSAligner = AlignerT<int8_fast, false, true>;
  using MSWordAligner = AlignerT<int16_fast, false, true>;
  using MSAlignerETE = AlignerT<int8_fast, true, true>;
  using MSWordAlignerETE = AlignerT<int16_fast, true, true>;


}

TEST_SUITE("Aligners");

TEST_CASE("Alignment") {

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
        reads.push_back("CCTT");
        reads.push_back("GGTT");
        reads.push_back("AAGG");
        reads.push_back("AACC");
        reads.push_back("AGGGT");
        reads.push_back("GG");
        reads.push_back("AAATTTA");
        reads.push_back("AAAGCCC");

        vargas::Results aligns;
        {
            vargas::Aligner a(7);
            aligns = a.align(reads, g.begin(), g.end());
        }
        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);
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

        // hisat like params
        vargas::Aligner a(10, 2, 6, 5, 3);
        vargas::Results aligns = a.align(reads, g.begin(), g.end());

        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);

        CHECK(aligns.max_score[8] == 12);
        CHECK(aligns.max_pos[8] == 10);

        CHECK(aligns.max_score[9] == 8);
        CHECK(aligns.max_pos[9] == 10);
    }

    SUBCASE("Quality") {
        std::vector<std::string> reads = {"GGTCTA", "GGTCTA", "GGTCTA"};
        std::vector<std::vector<char>> quals = {{40,40,40,0,40,40},
                                                {40,40,40,10,40,40},
                                                {40,40,40,20,40,40}};
        vargas::ScoreProfile prof(2, 2, 10, 10);
        prof.mismatch_min = 2;
        prof.mismatch_max = 6;
        vargas::Results aligns;
        vargas::Aligner a(6, prof);
        a.align_into(reads, quals, g.begin(), g.end(), aligns, true);
        REQUIRE(aligns.size() == 3);
        CHECK(prof.penalty(0) == 2);
        CHECK(prof.penalty(10) == 3);
        CHECK(prof.penalty(20) == 4);
        CHECK(prof.penalty(30) == 5);
        CHECK(prof.penalty(40) == 6);
        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_score[1] == 7);
        CHECK(aligns.max_score[2] == 6);
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

        vargas::WordAligner a(7);
        vargas::Results aligns = a.align(reads, g.begin(), g.end());
        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);
    }

    SUBCASE("Scoring Scheme- Word") {

        /**
*     GGG
*    /   \
* AAA     TTTA
*    \   /
*     CCC(ref)
*/

        std::vector<std::string> reads;
        reads.push_back("CCTT");
        reads.push_back("GGTT");
        reads.push_back("AAGG");
        reads.push_back("AACC");
        reads.push_back("AGGGT");
        reads.push_back("GG");
        reads.push_back("AAATTTA");
        reads.push_back("AAAGCCC");
        reads.push_back("AAAGAGTTTA");
        reads.push_back("AAAGAATTTA");

        // hisat like params
        vargas::WordAligner a(10, 2, 6, 5, 3);
        vargas::Results aligns = a.align(reads, g.begin(), g.end());

        CHECK(aligns.max_score[0] == 8);
        CHECK(aligns.max_pos[0] == 8);

        CHECK(aligns.max_score[1] == 8);
        CHECK(aligns.max_pos[1] == 8);

        CHECK(aligns.max_score[2] == 8);
        CHECK(aligns.max_pos[2] == 5);

        CHECK(aligns.max_score[3] == 8);
        CHECK(aligns.max_pos[3] == 5);

        CHECK(aligns.max_score[4] == 10);
        CHECK(aligns.max_pos[4] == 7);

        CHECK(aligns.max_score[5] == 4);
        CHECK(aligns.max_pos[5] == 6);

        CHECK(aligns.max_score[6] == 8);
        CHECK(aligns.max_pos[6] == 10);

        CHECK(aligns.max_score[7] == 8);
        CHECK(aligns.max_pos[7] == 4);

        CHECK(aligns.max_score[8] == 12);
        CHECK(aligns.max_pos[8] == 10);

        CHECK(aligns.max_score[9] == 8);
        CHECK(aligns.max_pos[9] == 10);
    }
}

TEST_CASE("Reverse strand") {
    vargas::Graph g;
    {
        vargas::Graph::Node n;
        n.set_as_ref();
        n.set_seq("ACGCGATCGACGATCGAACGATCGATGCCAGTGC");
        n.set_endpos(33);
        g.add_node(n);
    }
    std::vector<std::string> reads = {
    "GCCAGTGC", // mp: 34, fwd
    "GCACTGGC" // mp: 34, rev
    };
    vargas::AlignerETE a(8);
    SUBCASE("Alignment") {
        auto res = a.align(reads, g.begin(), g.end(), false);
        REQUIRE(res.size() == 2);
        CHECK(res.max_pos[0] == 34);
        CHECK(res.max_pos[1] == 34);
        CHECK(res.max_strand[0] == vargas::Strand::FWD);
        CHECK(res.max_strand[1] == vargas::Strand::REV);
    }
}

TEST_CASE("Indels") {
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

TEST_CASE("End to End alignment") {
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

TEST_CASE("Target score") {
    vargas::Graph g;
    vargas::Graph::Node n;
    n.set_seq("AAAACCCCCCCCCCCCAAA"); // Length 19
    n.set_endpos(18);
    g.add_node(n);

    const std::vector<std::string> reads = {"AAAA"};
    vargas::Aligner aligner(4);
    auto res = aligner.align(reads, g.begin(), g.end());
    REQUIRE(res.size() == 1);
    CHECK(res.max_score[0] == 8);
    CHECK(res.sub_score[0] == 6);
    CHECK(res.max_pos[0] == 4);
    CHECK(res.sub_pos[0] == 19);
}

TEST_SUITE_END();

#endif //VARGAS_ALIGNMENT_H
