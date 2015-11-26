/**
 * Ravi Gaddipati
 * November 24, 2015
 * rgaddip1@jhu.edu
 *
 * ReadSim simulates reads given a gssw_graph or a Graph object.
 * Reads can either be from a single individual or a random walk
 * with adjustable error rates.
 *
 * readsim.h
 */

#ifndef VARGAS_READSIM_H
#define VARGAS_READSIM_H

#include "graph.h"
#include "readsource.h"
#include <regex>
#include <map>
#include <random>

namespace vargas {

struct SimParams {
  time_t seed = time(NULL);
  bool randWalk = false;
  double muterr = 0;
  double indelerr = 0;
  int32_t readLen = 100;
  int maxreads = 10000;
};

inline std::ostream &operator<<(std::ostream &os, const SimParams &p) {
  std::stringstream ss;
  ss << "Seed: " << p.seed
      << ", Mut err: " << p.muterr
      << ", Indel Err: " << p.indelerr
      << ", read Len: " << p.readLen
      << ", max reads: " << p.maxreads
      << ", random walk? " << p.randWalk;
  os << ss.str();
  return os;
}

class ReadSim: public ReadSource {
 public:
  ReadSim() {
    srand(p.seed);
  }
  ReadSim(Graph &g) {
    ReadSim();
    setGraph(g);
  }
  ReadSim(SimParams param) {
    setParams(param);
  }
  ~ReadSim() {
    for (auto e : logs) {
      (*(e.second)).close();
      delete (e.second);
    }
  }
  void setGraph(Graph &g) {
    graph = g.getGSSWGraph();
  }
  void setGraph(gssw_graph *g) {
    graph = g;
  }

  Read &getRead() { return read; }
  bool updateRead();

  // Set the mutation error
  void setMuterr(double err) { p.muterr = err; }

  // Set the indel error
  void setIndelerr(double err) { p.indelerr = err; }

  // Set the length of the read
  void setReadLen(int32_t len) { p.readLen = len; }

  // Don't simulate from a specific individual
  void setRandWalk(bool randwalk) { p.randWalk = randwalk; }

  // Number of reads to generate for each regex
  void setNumReads(uint32_t num) { p.maxreads = num; }

  void setParams(SimParams param) {p = param;}

  // Add a regex, generates maxreads of each
  void addRegex(std::string regex, std::string file) {
    regexps.push_back(regex);
    logs.emplace(regex.c_str(), new std::ofstream(file));
    counters.emplace(regex.c_str(), 0);
    if (!logs[regex]->good()) throw std::invalid_argument("Error opening file: " + file);
    *(logs[regex]) << "#" << regex << std::endl
        << '#' << p << std::endl;
  }
  void clearRegexps() { regexps.clear(); }

 protected:
  SimParams p;
  std::vector<std::string> regexps;
  std::map<std::string, int> counters;
  std::map<std::string, std::ofstream *> logs;
  gssw_graph *graph = NULL;
  int totalreads = 0;

  void generateRead();

};

}


#endif //VARGAS_READSIM_H
