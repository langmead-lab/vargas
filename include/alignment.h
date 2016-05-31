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
#include "simdpp/simd.h"


namespace vargas {

/**
 * Stores a Read and associated alignment information. Positions are
 * 1 indexed.
 * @param read Read that was aligned
 * @param opt_score Best alignment score
 * @param opt_align_end Best alignment position, indexed w.r.t the last base of read
 * @param opt_count Number of alignments that tied for the best score
 * @param sub_score Second best alignment score
 * @param sub_align_end Second best alignment position
 * @param sub_count Number of second-best alignments
 * @param corflag 0 if alignment matches read origin, 1 of second best, 2 otherwise
 */
struct Alignment {
  Read read;

  uint16_t opt_score;
  int32_t opt_align_end;
  int32_t opt_count;

  uint16_t sub_score;
  int32_t sub_align_end;
  int32_t sub_count;

  int8_t corflag;


  Alignment()
      : opt_score(0), opt_align_end(-1), opt_count(-1), sub_score(0), sub_align_end(-1), sub_count(-1),
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

      this->read.end_pos = std::stoi(splitLine[0]);
      this->read.indiv = std::stoi(splitLine[1]);
      this->read.sub_err = std::stoi(splitLine[2]);
      this->read.var_nodes = std::stoi(splitLine[3]);
      this->read.var_bases = std::stoi(splitLine[4]);

      this->opt_score = (uint16_t) std::stoi(splitLine[5]);
      this->opt_align_end = std::stoi(splitLine[6]);
      this->opt_count = std::stoi(splitLine[7]);
      this->sub_score = (uint16_t) std::stoi(splitLine[8]);
      this->sub_align_end = std::stoi(splitLine[9]);
      this->sub_count = std::stoi(splitLine[10]);
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
  os << a.read << ',' << a.opt_score << ',' << a.opt_align_end << ',' << a.opt_count
      << ',' << a.sub_score << ',' << a.sub_align_end << ',' << a.sub_count
      << ',' << int32_t(a.corflag);
  return os;
}

template<int block_size>
class BatchAligner {
 public:
  BatchAligner() { }

  //private:
  simdpp::uint8<block_size> _packaged_reads = NULL;

  void _package_reads(const std::vector<Read> &reads) {
    // TODO optimize
  }

};
TEST_CASE ("Batch Aligner") {
  simdpp::aligned_allocator<simdpp::uint8<4>, 4> a;
  simdpp::uint8<4> *_packaged_reads = a.allocate(4);
}

}

#endif //VARGAS_ALIGNMENT_H
