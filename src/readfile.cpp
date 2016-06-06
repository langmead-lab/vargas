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
#include "utils.h"
#include "readfile.h"

bool vargas::ReadFile::update_read() {
  if (!_read_file) throw std::invalid_argument("No _read_file specified.");
  if (!std::getline(_read_file, line)) return false;
  if (line.at(0) == '#') {
    header += "\n" + line;
    return update_read();
  }

  unsigned long delim = line.find('#');
  if (delim == std::string::npos) {
    read.read = line;
    read.end_pos = 0;
    read.indiv = -1;
    read.sub_err = -1;
    read.indel_err = -1;
    read.var_nodes = -1;
    read.var_bases = -1;
    read.read_num = seq_to_num(line);
  } else {
    read.read = line.substr(0, delim);
    split(line.substr(delim + 1), ',', _split_meta);
    if (_split_meta.size() != 6) {
      std::cerr << "Unexpected number of fields." << std::endl;
      read.end_pos = 0;
      read.indiv = -1;
      read.sub_err = -1;
      read.indel_err = -1;
      read.var_nodes = -1;
      read.var_bases = -1;
      return false;
    }

    read.end_pos = uint32_t(std::stoi(_split_meta[0]));
    read.indiv = std::stoi(_split_meta[1]);
    read.sub_err = std::stoi(_split_meta[2]);
    read.indel_err = std::stoi(_split_meta[3]);
    read.var_nodes = std::stoi(_split_meta[4]);
    read.var_bases = std::stoi(_split_meta[5]);
    read.read_num = seq_to_num(read.read);
  }
  return true;
}

void vargas::ReadFile::resume_from(std::string read) {
  do {
    if (!update_read()) {
      std::cerr << "Warning: read not found in reads file. Starting from beginning." << std::endl;
      _read_file.clear();
      _read_file.seekg(0, std::ios::beg);
      break;
    }
  } while (get_read().read != read);
}