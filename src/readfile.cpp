/**
 * Ravi Gaddipati
 * November 24, 2015
 * rgaddip1@jhu.edu
 *
 * Wrapper for a reads file that loads reads and meta information
 * into a Read structure.
 *
 * readfile.cpp
 */

#include <string>
#include <fstream>
#include <stdexcept>
#include <vector>
#include "../include/utils.h"
#include "../include/readfile.h"

bool vargas::ReadFile::updateRead() {
  if (!readfile) throw std::invalid_argument("No readfile specified.");
  if (!std::getline(readfile, line)) return false;
  if (line.at(0) == '#') {
    header += "\n" + line;
    return updateRead();
  }

  unsigned long delim = line.find('#');
  if (delim == std::string::npos) {
    read.read = line;
    read.readEnd = 0;
    read.indiv = -1;
    read.numSubErr = -1;
    read.numIndelErr = -1;
    read.numVarNodes = -1;
    read.numVarBases = -1;
  } else {
    read.read = line.substr(0, delim);
    split(line.substr(delim + 1), ',', splitMeta);
    if (splitMeta.size() != 6) {
      std::cerr << "Unexpected number of fields." << std::endl;
      read.readEnd = 0;
      read.indiv = -1;
      read.numSubErr = -1;
      read.numIndelErr = -1;
      read.numVarNodes = -1;
      read.numVarBases = -1;
      return false;
    }

    read.readEnd = uint32_t(std::stoi(splitMeta[0]));
    read.indiv = std::stoi(splitMeta[1]);
    read.numSubErr = std::stoi(splitMeta[2]);
    read.numIndelErr = std::stoi(splitMeta[3]);
    read.numVarNodes = std::stoi(splitMeta[4]);
    read.numVarBases = std::stoi(splitMeta[5]);
  }
  return true;
}