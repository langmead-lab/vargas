/**
 * Ravi Gaddipati
 * November 24, 2015
 * rgaddip1@jhu.edu
 *
 * Abstract class for objects that can be used as a
 * source of reads (i.e. ReadSim and ReadFile).
 *
 * readsource.h
 */

#ifndef VMATCH_READS_H
#define VMATCH_READS_H

#include <string>
#include <sstream>

namespace vargas {

struct Read {
  std::string read;
  uint32_t readEnd;
  int32_t indiv;
  int32_t numSubErr;
  int32_t numVarNodes;
  int32_t numVarBases;

};

inline std::ostream &operator<<(std::ostream &os, const Read &r) {
  std::stringstream ss;
  ss << r.read << '#' << r.readEnd << ',' << r.indiv
      << ',' << r.numSubErr << ',' << r.numVarNodes
      << ',' << r.numVarBases;
  os << ss.str();
  return os;
}

class ReadSource {

 public:

  ReadSource() { }
  virtual ~ReadSource() { }

  // Updates the read and returns the string representation
  virtual std::string get() {
    if (!updateRead()) {
      read.read = "";
    }
    return str();
  }

  // Returns a string representation
  virtual std::string str() {
    std::stringstream ss;
    Read r = getRead();
    ss << r.read << '#' << r.readEnd << ',' << r.indiv
        << ',' << r.numSubErr << ',' << r.numVarNodes
        << ',' << r.numVarBases;
    return ss.str();
  };

  // Get the current read object
  virtual Read &getRead() = 0;

  // Update the current read, return false if none are available
  virtual bool updateRead() = 0;

  std::ostream &operator<<(std::ostream &os) {
    os << get();
    return os;
  }

 protected:
  Read read;

};

}

#endif //VMATCH_READS_H
