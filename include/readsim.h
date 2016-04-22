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

/*
 * Provides a 'target' read profile.
 * @param indiv Individual the read was taken from.
 * @param numSubErr Target number of substitiution errors introduced.
 * @param numVarNodes Target number of variant nodes the read traverses.
 * @param numVarBases Target number of bases that are in variant nodes.
 * @param numIndelErr Target number of insertions and deletions introduced.
 * @param readLen Target read length.
 */
struct ReadProfile {
  int32_t numSubErr = -1;
  int32_t numVarNodes = -1;
  int32_t numVarBases = -1;
  int32_t numIndelErr = -1;
  int32_t readLen = -1;

  ReadProfile() { }

  ReadProfile(const ReadProfile &p) {
    numSubErr = p.numSubErr;
    numVarNodes = p.numVarNodes;
    numVarBases = p.numVarBases;
    numIndelErr = p.numIndelErr;
    readLen = p.readLen;
  }

  /**
   * Checks if a Read matches the profile.
   * @return True if it matches.
   */
  bool matches(Read r) {
    if (numSubErr >= 0 && r.numSubErr != numSubErr) return false;
    if (numVarNodes >= 0 && r.numVarNodes != numVarNodes) return false;
    if (numVarBases >= 0 && r.numVarBases != numVarBases) return false;
    if (numIndelErr >= 0 && r.numIndelErr != numIndelErr) return false;
    if (readLen >= 0 && r.read.length() != readLen) return false;
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
  os << "numSubErr=" << p.numSubErr
      << ",numIndelErr=" << p.numIndelErr
      << ",numVarNodes=" << p.numVarNodes
      << ",numVarBases=" << p.numVarBases
      << ",readLen=" << p.readLen;
  return os;
}

/**
 * Parameters for the read simulator.
 * @param seed Use a predefined random seed
 * @param randWalk Ignore individuals, choose a random path
 * @param muterr Mutation error rate
 * @param indelerr Insertion and Ddletion error rate.
 * @param readLen Target read length. May be shorter if reach graph end.
 * @param ambiguity Limit on how many ambigous bases we allow.
 */
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
    setParams(param);
  }
  ~ReadSim() {
    for (auto &e : logs) {
      (*(e.second)).close();
      delete (e.second);
      delete (e.first);
    }
  }

  /**
   * Generates reads from the provided graph.
   * @param g Graph to pull reads from.
   */
  void setGraph(Graph &g) {
    graph = g.getGSSWGraph();
  }
  void setGraph(gssw_graph *g) {
    graph = g;
  }

  /**
   * Runs the sim until all provided profiles are fufilled.
   */
  void populateProfiles();

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

  void setParams(SimParams param) { p = param; }

  /**
   * Add a profile of read types we want.
   * @param prof Target read profile.
   * @param file Output filename.
   */
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
  void clearProfiles() { readProfiles.clear(); }


 protected:

  SimParams p;
  std::vector<ReadProfile *> readProfiles;
  std::map<ReadProfile *, int> counters;
  std::map<ReadProfile *, std::ofstream *> logs;
  gssw_graph *graph = NULL;
  int totalreads = 0;

};

}


#endif //VARGAS_READSIM_H
