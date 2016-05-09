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

/**
 * Stores a Read and associated alignment information. Positions are
 * 1 indexed.
 * @param read Read that was aligned
 * @param optScore Best alignment score
 * @param optAlignEnd Best alignment position, indexed w.r.t the last base of read
 * @param optCount Number of alignments that tied for the best score
 * @param subOptScore Second best alignment score
 * @param subOptAlignEnd Second best alignment position
 * @param subOptCount Number of second-best alignments
 * @param corflag 0 if alignment matches read origin, 1 of second best, 2 otherwise
 */
struct Alignment {
  Read read;

  uint16_t optScore;
  int32_t optAlignEnd;
  int32_t optCount;

  uint16_t subOptScore;
  int32_t subOptAlignEnd;
  int32_t subOptCount;

  int8_t corflag;


  Alignment()
      : optScore(0), optAlignEnd(-1), optCount(-1), subOptScore(0), subOptAlignEnd(-1), subOptCount(-1),
        corflag(-1) { }

  /**
   * Create an alignment with a read and associated meta information.
   * If meta information does not match the expected format, return after
   * populating the raw read.
   * @param line Read
   */
  Alignment(std::string line) {
    Alignment();
    std::vector<std::string> splitLine = split(line, '#');
    // We have the read string
    if (splitLine.size() > 0) {
      this->read.read = splitLine[0];
    }

    // We have meta info
    if (splitLine.size() > 1) {
      splitLine = split(splitLine[1], ',');
      if (splitLine.size() != 12) {
        // Unexpected format
        return;
      }

      this->read.readEndPos = std::stoi(splitLine[0]);
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
  }

};

/**
 * Print the alignment to os. This ordering matches the way the alignment is parsed
 * from a string.
 * @param os Output stream
 * @param an Alignment output
 */
inline std::ostream &operator<<(std::ostream &os, const Alignment &a) {
  os << a.read << ',' << a.optScore << ',' << a.optAlignEnd << ',' << a.optCount
      << ',' << a.subOptScore << ',' << a.subOptAlignEnd << ',' << a.subOptCount
      << ',' << int32_t(a.corflag);
  return os;
}

}

#endif //VARGAS_ALIGNMENT_H
