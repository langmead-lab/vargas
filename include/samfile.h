/**
 * Ravi Gaddipati
 * April 22, 2016
 * rgaddip1@jhu.edu
 *
 * Interface to work with SAM files.
 *
 * samfile.h
 */

#ifndef VARGAS_SAMFILE_H
#define VARGAS_SAMFILE_H

#include "readsource.h"
#include <string>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <map>

namespace vargas {

class SAMFile: public ReadSource {

 public:
  /**
   * Header reference sequence line, @SQ
   * @param name SN, reference seq name
   * @param length LN, ref seq length
   * @param genAsm AS, Genome Assembly identifier
   * @param md5 M5, checksum of seq
   * @param species SP, species
   * @param ur UR, URI of sequence, either web or local addr
   */
  struct ReferenceSequence {
    std::string raw;
    std::string name;
    int length;
    std::string genAsm;
    std::string md5;
    std::string species;
    std::string ur;
  };


  /**
   * Header read group line, @RG
   * @param id ID, Read group identifier
   * @param name CN, Name of sequencing center
   * @param desc DS
   * @param date DT, date run was produced
   * @param flowOrder FO
   * @param keySeq KS, array of base that corr. to key seq
   * @param lib LB, library
   * @param prog PG, Program used to process
   * @param median PI, predicted median insert size
   * @param platform PL, Platform/tech used
   * @param model PM, Platform model
   * @param unit PU, Platform unit
   * @param sample SM, sample or pool name
   */
  struct ReadGroup {
    std::string raw;
    std::string id;
    std::string name;
    std::string desc;
    std::string date;
    std::string flowOrder;
    std::string keySeq;
    std::string lib;
    std::string prog;
    std::string median;
    std::string platform;
    std::string model;
    std::string unit;
    std::string sample;
  };

  /**
   * Program header line, PG
   * @param id ID, program record identifier
   * @param name PN, program name
   * @param commandLine CL, exec command line
   * @param prevPrg PP, previous program applied
   * @param desc DS, description
   * @param version VN, program version
   */
  struct Program {
    std::string raw;
    std::string id;
    std::string name;
    std::string commandLine;
    std::string prevPrg;
    std::string desc;
    std::string version;
  };

  struct Flag {
    bool multipleSegments = false;
    bool properlyAligned = false;
    bool unmapped = false;
    bool nextSegmentUnmapped = false;
    bool reverseComplemented = false;
    bool nextReverseComplemented = false;
    bool firstSeg = false;
    bool lastSeg = false;
    bool secondaryAlignment = false;
    bool notPassingFilters = false;
    bool duplicate = false;
    bool supplementaryAlignment = false;

    Flag() { }
    Flag(unsigned int f) {
      multipleSegments = f & 0x1;
      properlyAligned = f & 0x2;
      unmapped = f & 0x4;
      nextSegmentUnmapped = f & 0x8;
      reverseComplemented = f & 0x10;
      nextReverseComplemented = f & 0x20;
      firstSeg = f & 0x40;
      lastSeg = f & 0x80;
      secondaryAlignment = f & 0x100;
      notPassingFilters = f & 0x200;
      duplicate = f & 0x400;
      supplementaryAlignment = f & 0x800;
    }
  };

  /**
   * Alignment record.
   * @param qname Query template name
   * @param flag bitwise flag
   * @param rname reference seq name
   * @param pos leftmost mapping position
   * @param mapq Mapping quality
   * @param cigar cigar string
   * @param rnext ref name of mate/next read
   * @param pnext position of mate/next read
   * @param tlen template length
   * @param seq segment sequence
   * @param qual ASCII Phred-scaled quality+33
   */
  struct Alignment {
    std::string raw;
    std::string qname;
    Flag flags;
    std::string rname;
    unsigned long pos;
    int mapq;
    std::string cigar;
    std::string rnext;
    unsigned long pnext;
    int tlen;
    std::string seq;
    std::string qual;

    std::map<std::string, std::string> optionalFields;
  };


  SAMFile() { }

  SAMFile(std::string file) {
    open(file);
  }

  ~SAMFile() {
    if (input.is_open()) input.close();
  }

  void open(std::string file) {
    if (input.is_open()) input.close();
    input.open(file);
    if (!input.good()) {
      throw std::invalid_argument("Invalid file: " + file);
    }
    init();
  }

  Read &get_read() override;
  std::string get_header() const override;
  bool update_read() override;

 protected:
  // Load the SAM file
  void init();
  void parseHeaderLine(std::string &line);
  Alignment parseAlignment(std::string &line);
  std::string reverseComplement(std::string &seq);


  // Data structures to hold the header
  std::vector<std::string> _comments; // @CO lines
  std::map<std::string, ReferenceSequence> refSeqs; // Maps a sequence name to its record
  std::unordered_map<std::string, ReadGroup> readGroups; // Maps a RG ID to a ReadGroup
  std::unordered_map<std::string, Program> programs; // Maps a program ID to a program

  Alignment currAlignment;
  Read currRead;
  std::ifstream input;
};

inline std::ostream &operator<<(std::ostream &os, const SAMFile::ReferenceSequence &rsq) {
  os << rsq.raw;
  return os;
}
inline std::ostream &operator<<(std::ostream &os, const SAMFile::Program &r) {
  os << r.raw;
  return os;
}
inline std::ostream &operator<<(std::ostream &os, const SAMFile::Alignment &r) {
  os << r.raw;
  return os;
}
inline std::ostream &operator<<(std::ostream &os, const SAMFile::ReadGroup &r) {
  os << r.raw;
  return os;
}

}


#endif //VARGAS_SAMFILE_H
