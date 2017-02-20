/**
 * @author Ravi Gaddipati
 * @date Feb 12, 2017
 * rgaddip1@jhu.edu
 *
 * @brief
 * Computes a traceback given a graph and location of the max score.
 * Smith-Waterman is computed from the max score position and extended
 * backwards.
 *
 * @details
 * Function max_extension_width determines how far back the traceback
 * will be computed before failing.
 */

#ifndef VARGAS_TRACEBACK_H
#define VARGAS_TRACEBACK_H

#include "graph.h"
#include "scoring.h"
#include "doctest.h"

#include <vector>
#include <algorithm>

namespace vargas {

  constexpr unsigned max_extension_width(const unsigned read_len) {return read_len * 2;}

  /**
   * @brief
   * Column major Smith-Waterman grid. SWGrid[0] is column zero and should
   * be Seed type.
   * @tparam CellT
   */
  template<typename CellT>
  class SWGrid {
    public:

      struct Cell {
          /**
           * Operation relative to the read.
           * Ir: Gap in the read
           * Iq: Gap in the reference
           * M: Match or mismatch
           */
          Cell() = default;
          virtual ~Cell() {}

          CellT scores[3]; // H, Insertion into query, Insertion into reference
          char origin; // Which matrix arrived from, 0 for H, 1 for Iq, 2 for Ir
      };

      /**
       * @brief
       * A seed is the same as a Cell, but also needs to keep track of which node it came from.
       */
      struct Seed : public Cell {
          unsigned id[3];
      };

      SWGrid() = default;
      SWGrid(unsigned height, unsigned width) {
          resize(height, width);
      }

      /**
       * @details
       * Padded with +1 in both dimensions. Column 0 should be seed.
       * @param height read length
       * @param width
       */
      void resize(unsigned height, unsigned width) {
          _grid.resize(width);
          std::for_each(_grid.begin(), _grid.end(), [height](std::vector<Cell> &v){
              v.resize(height + 1);
              v[0] = std::numeric_limits<CellT>::min();
          });
          _seeds.resize(height + 1);
      }

      /**
       * @param i column index
       * @return Column i
       */
      std::vector<Cell> &operator[](int i) {
          return _grid[i];
      }

      Seed &seed(unsigned i) {
          return _seeds[i];
      }

      Seed &seed(unsigned i) const {
          return _seeds.at(i);
      }

      std::vector<std::vector<Cell>>::iterator begin() {
          return _grid.begin();
      }

      std::vector<std::vector<Cell>>::iterator end() {
          return _grid.end();
      }

      std::vector<Cell> &back() {
          return _grid.back();
      }

      std::vector<Cell> &back() const {
          return _grid.back();
      }

    private:
      std::vector<std::vector<Cell>> _grid; // Fill col first, so grid[0] is col 0
      std::vector<Seed> _seeds;
  };

  /**
   * @brief
   * Computes a full traceback given a graph and the location of the maximal scoring cell.
   * @details
   * The Smith-Waterman grid is filled in from the end of the read and the maximal position,
   * and traced back.
   * @tparam CellT Score cell type
   * @tparam END_TO_END
   */
  template<typename CellT, bool END_TO_END>
  class Tracer {
    public:
      Tracer(const unsigned read_len, const ScoreProfile &prof)
      : _prof(prof), _read_len(read_len), _max_w(max_extension_width(_read_len)) {}

      std::string traceback(std::string read, const Graph &graph, unsigned maxs_pos) {
          auto curr_n = std::find_if(graph.rbegin(), graph.rend(),
                                     [maxs_pos](const Graph::Node &n){return n.begin_pos() < maxs_pos;});
          std::reverse(read.begin(), read.end());

          do {
              auto seq_beg = curr_n->rbegin(), seq_end;
              if (curr_n->end_pos() > maxs_pos) {
                  // If node extends past max score position
                  seq_beg += curr_n->end_pos() - maxs_pos;
              }
              if (curr_n->begin_pos() < maxs_pos - _max_w) {
                  // If node starts before max score position minus max extension width
                  seq_end = curr_n->rbegin() + curr_n->end_pos() - maxs_pos - _max_w;
              }
              else seq_end = curr_n->rend();
              auto &ref = (_id_to_sw[curr_n->id()] = SWGrid(_read_len, std::distance(seq_beg, seq_end)));

              _seed_grid(ref, curr_n.outgoing());
              _fill_grid(ref, read, seq_beg, seq_end);

              ++curr_n;

          } while(curr_n->begin_pos() > maxs_pos + _max_w);
      }

    private:
      std::unordered_map<unsigned, SWGrid> _id_to_sw;
      ScoreProfile _prof;
      const unsigned _read_len, _max_w;

      void _seed_grid(SWGrid &grid, const std::vector<unsigned> &outgoing) {
          for (auto pid : outgoing) {
              if (!_id_to_sw.count(pid)) continue;
              for (unsigned i = 1; i <= _read_len; ++i) {
                  auto &s = grid.seed(i);
                  const SWGrid::Cell &n_cell = _id_to_sw[pid].back()[i];
                  if (n_cell.scores[0] > s.scores[0]) {
                      s.scores[0] = n_cell.scores[0];
                      s.id[0] = pid;
                  }
                  if (n_cell.scores[1] > s.scores[1]) {
                      s.scores[1] = n_cell.scores[1];
                      s.id[1] = pid;
                  }
                  if (n_cell.scores[2] > s.scores[2]) {
                      s.scores[2] = n_cell.scores[2];
                      s.id[2] = pid;
                  }
              }
          }
      }

      void _fill_grid(SWGrid &grid, const std::string &read,
                      const std::vector<rg::Base>::const_reverse_iterator &seq_begin,
                      const std::vector<rg::Base>::const_reverse_iterator &seq_end) {
          const unsigned len = std::distance(seq_begin, seq_end);

          CellT tmp;
          //TODO saturation

          // Col 0
          for (unsigned r = 1; r <= _read_len; ++r) {
              SWGrid::Cell &curr = grid[0][r];

              curr.scores[1] = grid.seed(r).scores[0] - _prof.read_gopen;
              tmp = grid.seed(r).scores[1] - _prof.read_gext;
              if (tmp > curr.scores[1]) curr.scores[1] = tmp;


              curr.scores[2] = grid[0][r - 1].scores[0] - _prof.ref_gopen;
              tmp = grid[0][r - 1].scores[2] - _prof.ref_gext;
              if (tmp > curr.scores[2]) curr.scores[2] = tmp;

              curr.scores[0] =


          }
          // Rest cols
          for (unsigned c = 1; c < len; ++c) {
              for (unsigned r = 1; r <= _read_len; ++r) {

              }
          }
      }

  };

}
#endif //VARGAS_TRACEBACK_H
