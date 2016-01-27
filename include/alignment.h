/**
 * Ravi Gaddipati
 * Jan 27, 2016
 * rgaddip1@jhu.edu
 *
 * Provides tools to interact with Alignments.
 *
 * alignment.h
 */

#ifndef VARGAS_ALIGNMENT_H
#define VARGAS_ALIGNMENT_H

#include "readsource.h"
#include "utils.h"


namespace vargas {

// struct to hold alignment results
struct Alignment {
  Read read;

  // Optimal alignment
  uint16_t optScore;
  int32_t optAlignEnd;
  int32_t optCount;

  // Suboptimal alignment
  uint16_t subOptScore;
  int32_t subOptAlignEnd;
  int32_t subOptCount;

  // alignment flag
  int8_t corflag;

  Alignment()
      : optScore(0), optAlignEnd(-1), optCount(-1), subOptScore(0), subOptAlignEnd(-1), subOptCount(-1),
        corflag(-1) { };
  Alignment(std::string line) {
    Alignment();
    std::vector<std::string> splitLine = split(line, '#');
    // We have the read string
    if (splitLine.size() > 0) {
      this->read.read = splitLine[0];
    }

    // We have meta info
    if (splitLine.size() > 1) {
      split(splitLine[1], ',', splitLine);
      if (splitLine.size() != 11) {
        // Unexpected format
        std::cerr << "Invalid alignment record." << std::endl;
        return;
      }

      this->read.readEnd = (uint32_t) std::stoi(splitLine[0]);
      this->read.indiv = std::stoi(splitLine[1]);
      this->read.numSubErr = std::stoi(splitLine[2]);
      this->read.numVarNodes = std::stoi(splitLine[3]);
      this->read.numVarBases = std::stoi(splitLine[4]);

      this->optScore = (uint16_t) std::stoi(splitLine[5]);
      this->optAlignEnd = std::stoi(splitLine[6]);
      this->optCount = std::stoi(splitLine[7]);
      this->subOptScore = (uint16_t) std::stoi(splitLine[8]);
      this->subOptAlignEnd = std::stoi(splitLine[9]);
      this->subOptCount = std::stoi(splitLine[10]);
      this->corflag = (int8_t) std::stoi(splitLine[11]);
    }
  };
};

}

#endif //VARGAS_ALIGNMENT_H
