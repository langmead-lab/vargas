//
// Created by gaddra on 11/24/15.
//

#ifndef VMATCH_READFILE_H
#define VMATCH_READFILE_H

#include <string>
#include <fstream>
#include <stdexcept>
#include "readsource.h"

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

  std::string getRead() { return read; }
  std::string getMeta() { return meta; }
  bool updateRead();

 protected:
  std::string line;
  std::ifstream readfile;

};

}

#endif //VMATCH_READFILE_H
