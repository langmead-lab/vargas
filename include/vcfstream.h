/**
 * Ravi Gaddipati
 * November 23, 2015
 * rgaddip1@jhu.edu
 *
 * vmatch::vcfstream is a filtering wrapper for a VCF file.
 * Variant records are parsed and filtered according to the ingroup.
 *
 * vcfstream.h
 */

#ifndef VMATCH_VCFSTREAM_H
#define VMATCH_VCFSTREAM_H

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include "../include/utils.h"
#include <vector>
#include <map>
#include <cstdlib>
#include <string>

namespace vmatch {

/** record of data that contains the parsed and filtered VCF **/
struct vcfrecord {
  ulong pos; // position of variant
  std::string ref; // ref is also included in indivs and freqs
  std::map<std::string, std::vector<uint32_t>> indivs; // Maps each alt to a list of indivs w/ that alt
  std::map<std::string, double_t> freqs; // Maps each alt to its allele frequency
};

// For printing a record
std::ostream &operator<<(std::ostream &os, const vcfrecord &vrec);

/**
 * Handles a VCF and parses the results. Filters variants based on the ingroup.
 */
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

  // get the next line in the VCF file and parse it
  bool getRecord(vcfrecord &vrecord);
  // Create a random ingroup.
  void createIngroup(int32_t percent = 100, long seed = NULL);
  // Define an ingroup
  void createIngroup(std::vector<uint32_t> &vec) {
    ingroup = vec;
    std::sort(ingroup.begin(), ingroup.end());
  };
  // Create a complement ingroup
  void createComplementIngroup(std::vector<uint32_t> vec);
  // Return a copy of the ingroup
  std::vector<uint32_t> getIngroup() const { return ingroup; };
  // Print a CSV of the ingroup individuals
  void printIngroup(std::ostream &os);
  // Update a record with the next filtered VCF line
  friend vmatch::vcfstream &operator>>(vcfstream &vstream, vcfrecord &vrec);

 protected:
  std::ifstream vcfFile; // source VCF file
  std::string currentRecord; // Most recent line read from VCF
  std::vector<std::string>
      splitTemp, // Used as a temp vector to split various items
      splitRecord, // The entire split currentRecord
      splitDiploid; // Used to split haplotypes
  std::vector<uint32_t> ingroup; // List of all individuals in the ingroup. Should always be sorted
  bool initilized = false;

  // Holds the locations of all of the fields in the VCF. Found during init.
  struct {
    int32_t chrom = -1;
    int32_t pos = -1;
    int32_t id = -1;
    int32_t ref = -1;
    int32_t alt = -1;
    int32_t qual = -1;
    int32_t filter = -1;
    int32_t info = -1;
    int32_t format = -1;
    int32_t indivOffset = -1; // indivOffset indicates the first pos i.e. indiv 0
    int32_t numIndivs = 0; // Each haplotype is an individual
  } fields;

 private:
  long seed = time(NULL); // ingroup generation seed
  std::vector<std::string> afSplit;
  std::vector<uint32_t> altIndivs;

  void initVCF();
  /**
   * Splits the currect record into a vector, including each diploid into haplotypes.
   * i.e. indivOffset = individual 0 haplotype 0
   *      indivOffset = individual 0 haplotype 1
  **/
  void splitCurrentRecord();
};

vmatch::vcfstream &operator>>(vcfstream &vstream, vcfrecord &vrec);

}

#endif //VMATCH_VCFSTREAM_H
