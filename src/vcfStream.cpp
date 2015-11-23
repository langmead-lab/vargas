//
// Created by gaddra on 11/22/15.
//

#include "../include/vcfStream.h"

void vmatch::vcfstream::initVCF() {
  // Skip the comments
  do { std::getline(vcfFile, currentRecord); } while (currentRecord.substr(0, 2) == "##");
  if (currentRecord.at(0) != '#') throw std::invalid_argument("Expected header beginning with #");

  // Convert to uppercase
  transform(currentRecord.begin(), currentRecord.end(), currentRecord.begin(), ::toupper);
  currentRecord = currentRecord.substr(1);
  std::vector<std::string> header = split(currentRecord, '\t');

  //Find locations
  fields.chrom = int32_t(std::find(header.begin(), header.end(), "CHROM") - header.begin());
  fields.pos = int32_t(std::find(header.begin(), header.end(), "POS") - header.begin());
  fields.id = int32_t(std::find(header.begin(), header.end(), "ID") - header.begin());
  fields.ref = int32_t(std::find(header.begin(), header.end(), "REF") - header.begin());
  fields.alt = int32_t(std::find(header.begin(), header.end(), "ALT") - header.begin());
  fields.qual = int32_t(std::find(header.begin(), header.end(), "QUAL") - header.begin());
  fields.filter = int32_t(std::find(header.begin(), header.end(), "FILTER") - header.begin());
  fields.info = int32_t(std::find(header.begin(), header.end(), "INFO") - header.begin());
  fields.format = int32_t(std::find(header.begin(), header.end(), "FORMAT") - header.begin());
  fields.indivOffset = fields.format + 1;
  fields.numIndivs = uint32_t((header.size() - fields.indivOffset) * 2);

  // Check required fields
  if (fields.pos < 0) throw std::invalid_argument("POS field not found.");
  if (fields.ref < 0) throw std::invalid_argument("ALT field not found.");
  if (fields.alt < 0) throw std::invalid_argument("REF field not found.");
  if (fields.info < 0) throw std::invalid_argument("INFO field not found.");
  if (fields.format < 0) throw std::invalid_argument("Currently only supports files with unphased genotype data.");

  initilized = true;

}

void vmatch::vcfstream::splitCurrentRecord() {
  splitTemp = split(currentRecord, '\t');
  splitRecord.clear();
  for (int i = 0; i < splitTemp.size(); ++i) {
    if (i >= fields.indivOffset) {
      splitHaplo = split(splitTemp[i], '|');
      splitRecord.push_back(splitHaplo[0]);
      splitRecord.push_back(splitHaplo[1]);
    } else {
      splitRecord.push_back(splitTemp[i]);
    }
  }
}


bool vmatch::vcfstream::getRecord(vmatch::vcfrecord &vrecord) {
  if (!initilized) {
    throw std::invalid_argument("VCF file not provided.");
  }
  if (!getline(vcfFile, currentRecord)) {
    return false;
  }
  splitCurrentRecord();
  vrecord.pos = ulong(atol(splitRecord[fields.pos].c_str()));
  vrecord.ref = splitRecord[fields.ref].c_str();
  vrecord.alts.clear();
  splitTemp = split(splitRecord[fields.alt], ',');
  std::vector<uint32_t> altIndivs;
  for (uint32_t i = 0; i < splitTemp.size(); i++) {
    altIndivs.clear();
    for (int32_t d = fields.indivOffset; d < fields.numIndivs + fields.indivOffset; ++d) {
      if (atoi(splitRecord[d].c_str()) == i + 1) {
        // Check if its in the ingroup
        if (std::find(ingroup.begin(), ingroup.end(), d) != ingroup.end() || ingroup.size() == 0)
          altIndivs.push_back(d);
      }
    }
    if (altIndivs.size() > 0)
      vrecord.alts.emplace(splitTemp[i].c_str(), altIndivs);
  }
  return true;
}