/**
 * Ravi Gaddipati
 * April 22, 2016
 * rgaddip1@jhu.edu
 *
 * Abstract class for objects that can be used as a
 * source of reads (i.e. ReadSim and ReadFile).
 *
 * readsource.h
 */

#ifndef VARGAS_READS_H
#define VARGAS_READS_H

#include <string>
#include <sstream>

namespace vargas {

/**
 * Struct to represent a Read.
 * @param read base sequence.
 * @param readEnd position of last base in seq.
 * @param indiv Individual the read was taken from.
 * @param numSubErr Number of substitiution errors introduced.
 * @param numVarNodes Number of variant nodes the read traverses.
 * @param numVarBases Number of bases that are in variant nodes.
 * @param numIndelErr Number of insertions and deletions introduced.
 */
struct Read {
  std::string read;
  int32_t readEnd;
  int32_t indiv;
  int32_t numSubErr;
  int32_t numVarNodes;
  int32_t numVarBases;
  int32_t numIndelErr;

};

inline std::ostream &operator<<(std::ostream &os, const Read &r) {
  std::stringstream ss;
  ss << r.read
      << '#' << r.readEnd
      << ',' << r.indiv
      << ',' << r.numSubErr
      << ',' << r.numIndelErr
      << ',' << r.numVarNodes
      << ',' << r.numVarBases;
  os << ss.str();
  return os;
}


class ReadSource {

 public:

  ReadSource() { }
  virtual ~ReadSource() { }

  // Updates the read and returns the string representation
  virtual std::string updateAndGet() {
    if (!updateRead()) {
      read.read = "";
    }
    return toString();
  }

  // Returns a string representation
  virtual std::string toString() {
    std::stringstream ss;
    Read r = getRead();
    ss << r.read
        << '#' << r.readEnd
        << ',' << r.indiv
        << ',' << r.numSubErr
        << ',' << r.numIndelErr
        << ',' << r.numVarNodes
        << ',' << r.numVarBases;
    return ss.str();
  };

  // Get the current read object
  virtual Read &getRead() = 0;

  // Get read file header
  virtual std::string getHeader() const = 0;

  // Update the current read, return false if none are available
  virtual bool updateRead() = 0;

  inline std::ostream &operator<<(std::ostream &os) {
    os << updateAndGet();
    return os;
  }


 protected:
  Read read;
  std::string header = "";

};

}

#endif //VARGAS_READS_H
