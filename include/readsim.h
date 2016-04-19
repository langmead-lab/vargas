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

#include "readsource.h"

namespace vargas {

struct ReadProfile {
  //TODO add numIndelErr to prof and Read
  int32_t indiv = -1;
  int32_t numSubErr = -1;
  int32_t numVarNodes = -1;
  int32_t numVarBases = -1;
    int32_t numIndelErr = -1;

  ReadProfile() { }

  ReadProfile(const ReadProfile &p) {
    indiv = p.indiv;
    numSubErr = p.numSubErr;
    numVarNodes = p.numVarNodes;
    numVarBases = p.numVarBases;
    numIndelErr = p.numIndelErr;
  }

  bool matches(Read r) {
    if (indiv >= 0 && r.indiv != indiv) return false;
    if (numSubErr >= 0 && r.numSubErr != numSubErr) return false;
    if (numVarNodes >= 0 && r.numVarNodes != numVarNodes) return false;
    if (numVarBases >= 0 && r.numVarBases != numVarBases) return false;
    if (numIndelErr >= 0 && r.numIndelErr != numIndelErr) return false;
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
      << ",numVarBases=" << p.numVarBases
      << ",numIndelErr=" << p.numIndelErr;
  return os;
}


struct SimParams {
  time_t seed = time(NULL);
  bool randWalk = false;
  double muterr = 0;
  double indelerr = 0;
  int32_t readLen = 100;
  int maxreads = 10000;
  int ambiguity = 0;

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
    srand(p.seed);
    setGraph(g);
  }
  ReadSim(SimParams param) {
    srand(p.seed);
    setParams(param);
  }
  ~ReadSim() {
    for (auto &e : logs) {
      (*(e.second)).close();
      delete (e.second);
      delete (e.first);
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
    ReadProfile *p = new ReadProfile(prof);
    std::ofstream *os = new std::ofstream(file);
    readProfiles.push_back(p);
    logs[p] = os;
    counters[p] = 0;
    if (!os->good()) throw std::invalid_argument("Error opening file: " + file);
    *(logs[p]) << "#" << (*p) << std::endl
        << '#' << this->p << std::endl;
  }
  void clearRegexps() { readProfiles.clear(); }


 protected:

  SimParams p;
  std::vector<ReadProfile *> readProfiles;
  std::map<ReadProfile *, int> counters;
  std::map<ReadProfile *, std::ofstream *> logs;
  gssw_graph *graph = NULL;
  int totalreads = 0;

  void generateRead();

};

}


#endif //VARGAS_READSIM_H
