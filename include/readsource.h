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
#include "utils.h"

namespace vargas {

/**
 * Struct to represent a Read.
 * @param read base sequence.
 * @param end_pos position of last base in seq.
 * @param indiv Individual the read was taken from.
 * @param sub_err Number of substitiution errors introduced.
 * @param var_nodes Number of variant nodes the read traverses.
 * @param var_bases Number of bases that are in variant nodes.
 * @param indel_err Number of insertions and deletions introduced.
 */
struct Read {
  std::string read;
  std::vector<uchar> read_num;
  int32_t end_pos;
  int32_t indiv;
  int32_t sub_err;
  int32_t var_nodes;
  int32_t var_bases;
  int32_t indel_err;

};

inline std::ostream &operator<<(std::ostream &os, const Read &r) {
  std::stringstream ss;
  ss << r.read
      << '#' << r.end_pos
      << ',' << r.indiv
      << ',' << r.sub_err
      << ',' << r.indel_err
      << ',' << r.var_nodes
      << ',' << r.var_bases;
  os << ss.str();
  return os;
}


class ReadSource {

 public:

  ReadSource() { }
  virtual ~ReadSource() { }

  /**
   * Updates the stored and and returns the read.
   */
  virtual std::string update_and_get() {
    if (!update_read()) {
      read.read = "";
    }
    return to_string();
  }

  // Returns a string representation
  virtual std::string to_string() {
    std::stringstream ss;
    Read r = get_read();
    ss << r.read
        << '#' << r.end_pos
        << ',' << r.indiv
        << ',' << r.sub_err
        << ',' << r.indel_err
        << ',' << r.var_nodes
        << ',' << r.var_bases;
    return ss.str();
  };

  // Get the current read object
  virtual Read &get_read() = 0;

  // Get read file header
  virtual std::string get_header() const = 0;

  // Update the current read, return false if none are available
  virtual bool update_read() = 0;

  const std::vector<Read> &get_batch(int size) {
    if (size <= 0) size = 1;
    if (_batch.size() != size) _batch.resize(size);
    for (int i = 0; i < size; ++i) {
      if (!update_read()) read = Read();
      _batch[i] = read;
    }
    return _batch;
  }

  inline std::ostream &operator<<(std::ostream &os) {
    os << update_and_get();
    return os;
  }


 protected:
  Read read;
  std::string header = "";
  std::vector<Read> _batch;

};

}

#endif //VARGAS_READS_H
