//
// Created by gaddra on 11/22/15.
//

#ifndef VMATCH_VCFSTREAM_H
#define VMATCH_VCFSTREAM_H

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include "../include/utils.h"
#include <vector>
#include <map>

namespace vmatch {

struct vcfrecord {
  vcfrecord() { }
  ulong pos;
  std::string ref;
  std::map<std::string, std::vector<uint32_t>> alts;
};

class vcfstream {
 public:
  vcfstream() { };
  vcfstream(std::string filename) {
    open(filename);
  }
  ~vcfstream() {
    vcfFile.close();
  }
  void open(std::string filename) {
    vcfFile.open(filename);
    if (!vcfFile.good()) throw std::invalid_argument("Invalid VCF File");
    initVCF();
  }
  bool getRecord(vcfrecord &vrecord);
 protected:
  std::ifstream vcfFile; // source VCF file
  std::string currentRecord;
  std::vector<std::string> splitTemp, splitRecord, splitHaplo;
  std::vector<uint32_t> ingroup;
  bool initilized = false;

  struct fieldPositions {
    int32_t chrom = -1;
    int32_t pos = -1;
    int32_t id = -1;
    int32_t ref = -1;
    int32_t alt = -1;
    int32_t qual = -1;
    int32_t filter = -1;
    int32_t info = -1;
    int32_t format = -1;
    int32_t indivOffset = -1;
    int32_t numIndivs = 0; // Each haplotype is an individual
  } fields;

  void initVCF();
  void splitCurrentRecord();
};
}

#endif //VMATCH_VCFSTREAM_H
