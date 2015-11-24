//
// Created by gaddra on 11/24/15.
//

#include "../include/readfile.h"

bool vmatch::ReadFile::updateRead() {
  if (!readfile) throw std::invalid_argument("No readfile specified.");
  if (!std::getline(readfile, line)) return false;
  if (line.at(0) == '#') return updateRead();
  unsigned long delim = line.find('#');
  if (delim == std::string::npos) {
    read = line;
    meta = "";
  } else {
    read = line.substr(0, delim);
    meta = line.substr(delim + 1);
  }
  return true;
}