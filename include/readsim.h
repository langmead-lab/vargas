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

#ifndef VMATCH_READSIM_H
#define VMATCH_READSIM_H

#include "graph.h"
#include "readsource.h"
#include <regex>
#include <map>
#include <random>

namespace vargas {

class ReadSim: public ReadSource {
 public:
  ReadSim() {
    srand(seed);
  }
  ReadSim(Graph &g) {
    ReadSim();
    setGraph(g);
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
  void setMuterr(double err) { muterr = err; }

  // Set the indel error
  void setIndelerr(double err) { indelerr = err; }

  // Set the length of the read
  void setReadLen(int32_t len) { readLen = len; }

  // Don't simulate from a specific individual
  void setRandWalk(bool randwalk) { randWalk = randwalk; }

  // Number of reads to generate for each regex
  void setNumReads(uint32_t num) { maxreads = num; }

  // Add a regex, generates maxreads of each
  void addRegex(std::string regex, std::string file) {
    regexps.push_back(regex);
    logs.emplace(regex.c_str(), new std::ofstream(file));
    counters.emplace(regex.c_str(), 0);
    *(logs[regex]) << "#" << regex << std::endl;
  }
  void clearRegexps() { regexps.clear(); }

 protected:
  time_t seed = time(NULL);
  std::vector<std::string> regexps;
  std::map<std::string, int> counters;
  std::map<std::string, std::ofstream *> logs;
  gssw_graph *graph = NULL;
  bool randWalk = false;
  double muterr = 0;
  double indelerr = 0;
  int32_t readLen = 100;
  int maxreads = 10000;
  int totalreads = 0;

  void generateRead();

};

}


#endif //VMATCH_READSIM_H
