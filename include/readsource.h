//
// Created by gaddra on 11/23/15.
//

#ifndef VMATCH_READS_H
#define VMATCH_READS_H

#include <string>
#include <sstream>

namespace vmatch {

class ReadSource {
 public:

  struct Read {
    std::string read;
    uint32_t readEnd;
    int32_t indiv;
    int32_t numSubErr;
    int32_t numVarNodes;
    int32_t numVarBases;
  };

  ReadSource() { }
  virtual ~ReadSource() { }

  virtual std::string get() {
    if (!updateRead()) {
      read.read = "";
    }
    return str();
  }
  virtual std::string str() {
    std::stringstream ss;
    Read r = getRead();
    ss << r.read << '#' << r.readEnd << ',' << r.indiv
        << ',' << r.numSubErr << ',' << r.numVarNodes
        << ',' << r.numVarBases;
    return ss.str();
  };
  virtual Read &getRead() = 0;
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
