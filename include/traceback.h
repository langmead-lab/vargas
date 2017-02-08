//
// Created by gaddra on 1/31/17.
//

#ifndef VARGAS_TRACEBACK_H
#define VARGAS_TRACEBACK_H

#include "graph.h"
#include "scoring.h"

#include <vector>
#include <algorithm>
#include <algorithm>

namespace vargas {

  template<typename CellT, bool END_TO_END>
  class Tracer {
    public:
      Tracer(const unsigned read_len, const ScoreProfile &prof)
      : _grid(read_len + (read_len / 2)), _prof(prof), _read_len(read_len) {
          std::for_each(_grid.begin(), _grid.end(), [](std::vector<CellT> &v) {
              v.resize(read_len + 1);
              v[0] = std::numeric_limits<CellT>::min();
          });
      }

      std::string traceback(const std::string &read, Graph &graph, unsigned maxs_pos) {
          auto subg = graph.subgraph(maxs_pos - _grid.size() + 1, maxs_pos);
          _seed seed;
          // We expect the entire alignment to be within this subgrpah, so seeding condition shouldn't matter
          std::fill(seed.S.begin(), seed.S.end(), std::numeric_limits<CellT>::min());
          std::fill(seed.I.begin(), seed.I.end(), std::numeric_limits<CellT>::min());

          for (auto n : subg) {

          }
      }


    private:
      enum class _Dir { M, I, D };

      struct _seed {
          _seed() = default;
          _seed(unsigned rl) : S(rl), I(rl), R(rl) {}
          std::vector<CellT> S, I;
          std::vector<_Dir> R;
      };

      void _get_seed(const std::vector<unsigned> &incoming, _seed &seed) {
          if (incoming.size() == 0) return;
          for (unsigned i = 1l i < _read_len + 1; ++i) {
              const auto &s = _seed_map.at(incoming[0]);
              seed.S[i] = s.S[i];
              seed.I[i] = s.I[i];
              for (unsigned p = 1; p < incoming.size(); ++p) {
                  seed.S[i] = std::max<CellT>(seed.S[i], s.S[i]);
                  seed.I[i] = std::max<CellT>(seed.S[i], s.I[i]);
              }
          }
      }


      std::vector<std::vector<CellT>> _grid; // cols x rows; _grid[0] is col 0
      std::vector<std::vector<_Dir>> _dir;
      std::unordered_map<unsigned, _seed> _seed_map;
      ScoreProfile _prof;
      const unsigned _read_len;
  };

}
#endif //VARGAS_TRACEBACK_H
