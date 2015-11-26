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
#include <map>
#include <random>

namespace vargas {

struct ReadProfile {
  //TODO add numIndelErr to prof and Read
  int32_t indiv = -1;
  int32_t numSubErr = -1;
  int32_t numVarNodes = -1;
  int32_t numVarBases = -1;

  bool matches(Read r) {
    if (indiv >= 0 && r.indiv != indiv) return false;
    if (numSubErr >= 0 && r.numSubErr != numSubErr) return false;
    if (numVarNodes >= 0 && r.numVarNodes != numVarNodes) return false;
    if (numVarBases >= 0 && r.numVarBases != numVarBases) return false;
    return true;
  }

};

inline bool operator==(Read &r, ReadProfile &p) {
  return p.matches(r);
}
inline bool operator!=(Read &r, ReadProfile &p) {
  return !p.matches(r);
}
inline bool operator==(ReadProfile &p, Read &r) {
  return p.matches(r);
}
inline bool operator!=(ReadProfile &p, Read &r) {
  return !p.matches(r);
}
inline std::ostream &operator<<(std::ostream &os, const ReadProfile &p) {
  os << "indiv=" << p.indiv
      << ",numSubErr=" << p.numSubErr
      << ",numVarNodes=" << p.numVarNodes
      << ",numVarBases=" << p.numVarBases;
  return os;
}

struct SimParams {
  time_t seed = time(NULL);
  bool randWalk = false;
  double muterr = 0;
  double indelerr = 0;
  int32_t readLen = 100;
  int maxreads = 10000;
};

inline std::ostream &operator<<(std::ostream &os, const SimParams &p) {
  os << "Seed: " << p.seed
      << ", Mut err: " << p.muterr
      << ", Indel Err: " << p.indelerr
      << ", read Len: " << p.readLen
      << ", max reads: " << p.maxreads
      << ", random walk? " << p.randWalk;
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
    ReadSim();
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
  std::string getHeader() const { return header; }

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
  void addProfile(ReadProfile &prof, std::string file) {
    readProfiles.push_back(prof);
    logs.emplace(&readProfiles.back(), new std::ofstream(file));
    counters.emplace(&readProfiles.back(), 0);
    if (!logs[&readProfiles.back()]->good()) throw std::invalid_argument("Error opening file: " + file);
    *(logs[&readProfiles.back()]) << "#" << readProfiles.back() << std::endl
        << '#' << p << std::endl;
  }
  void clearRegexps() { readProfiles.clear(); }

 protected:
  SimParams p;
  std::vector<ReadProfile> readProfiles;
  std::map<ReadProfile *, int> counters;
  std::map<ReadProfile *, std::ofstream *> logs;
  gssw_graph *graph = NULL;
  int totalreads = 0;

  void generateRead();

};

}


#endif //VARGAS_READSIM_H
