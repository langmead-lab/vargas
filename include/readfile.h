/**
 * Ravi Gaddipati
 * November 24, 2015
 * rgaddip1@jhu.edu
 *
 * Wrapper for a reads file that loads reads and meta information
 * into a Read structure.
 *
 * readfile.h
 */

#ifndef VMATCH_READFILE_H
#define VMATCH_READFILE_H

#include <string>
#include <fstream>
#include <stdexcept>
#include "readsource.h"
#include <vector>
#include "../include/utils.h"

namespace vmatch {

class ReadFile: public ReadSource {
 public:
  ReadFile() { }
  ReadFile(std::string file) {
    readfile.open(file);
    if (!readfile.good()) throw std::invalid_argument("Invalid read file.");
  }
  ~ReadFile() {
    if (readfile) readfile.close();
  }
  Read &getRead() { return read; }
  bool updateRead();

 protected:
  std::string line;
  std::ifstream readfile;
  std::vector<std::string> splitMeta;

};

}

#endif //VMATCH_READFILE_H
