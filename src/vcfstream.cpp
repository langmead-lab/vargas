//
// Created by gaddra on 11/22/15.
//

#include "../include/vcfstream.h"

void vmatch::vcfstream::initVCF() {
  // Skip the comments
  do { std::getline(vcfFile, currentRecord); } while (currentRecord.substr(0, 2) == "##");
  if (currentRecord.at(0) != '#') throw std::invalid_argument("Expected header beginning with #");

  // Convert to uppercase
  transform(currentRecord.begin(), currentRecord.end(), currentRecord.begin(), ::toupper);
  currentRecord = currentRecord.substr(1); // Remove leading #
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

  // Default ingroup is everyone
  createIngroup(100);

}

void vmatch::vcfstream::splitCurrentRecord() {
  split(currentRecord, '\t', splitTemp);
  splitRecord.clear();
  for (int i = 0; i < splitTemp.size(); ++i) {
    if (i >= fields.indivOffset) {
      split(splitTemp[i], '|', splitDiploid);
      splitRecord.push_back(splitDiploid[0]);
      splitRecord.push_back(splitDiploid[1]);
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
  vrecord.altIndivs.clear();
  vrecord.altAF.clear();
  split(splitRecord[fields.alt], ',', splitTemp);
  afSplit.clear();

  // Get the list of AF's
  split(splitRecord[fields.info], ';', afSplit);
  // Find and split the AF entry
  for (auto e : afSplit) {
    if (e.substr(0, 3) == "AF=") {
      split(e.substr(3), ',', afSplit);
      break;
    }
  }

  // Get the reference frequency
  double_t sumaltAF = 0;
  for (auto a : afSplit) {
    sumaltAF += atof(a.c_str());
  }
  vrecord.refFreq = 1 - sumaltAF;

  bool validAF = true;
  if (afSplit.size() != splitTemp.size()) {
    std::cerr << "Alternate and AF field lengths do not match at pos " << splitRecord[fields.pos] << ".\n";
    validAF = false;
  }

  // For each alternate allele
  for (int i = 0; i < splitTemp.size(); i++) {
    altIndivs.clear();
    // For each individual, check if it is in the ingroup
    for (int d = fields.indivOffset; d < fields.numIndivs + fields.indivOffset; ++d) {
      // a value of 0 indicates the reference is used, i+1 gives the alt number.
      if (atoi(splitRecord[d].c_str()) == i + 1) {
        // Check if its in the ingroup
        if (std::binary_search(ingroup.begin(), ingroup.end(), d)) altIndivs.push_back(d);
      }
    }
    if (altIndivs.size() > 0) {
      vrecord.altIndivs.emplace(splitTemp[i].c_str(), altIndivs);
      if (validAF)
        vrecord.altAF.emplace(splitTemp[i].c_str(), std::atof(afSplit[i].c_str()));
    }
  }
  return true;
}

void vmatch::vcfstream::createIngroup(int32_t percent, long seed) {
  if (!initilized) {
    throw std::invalid_argument("VCF file not provided.");
  }
  if (seed != NULL) this->seed = seed;
  ingroup.clear();

  if (percent == 100) {
    for (int i = fields.indivOffset; i < fields.numIndivs + fields.indivOffset; ++i) {
      ingroup.push_back(i);
      std::sort(ingroup.begin(), ingroup.end());
    }
  }
  else if (percent == 0) {
    return;
  } else {
    for (int i = fields.indivOffset; i < fields.numIndivs + fields.indivOffset; ++i) {
      if (rand() % 10000 < percent * 100) ingroup.push_back(i);
    }
    std::sort(ingroup.begin(), ingroup.end());
  }
}

void vmatch::vcfstream::createComplementIngroup(std::vector<uint32_t> vec) {
  std::sort(vec.begin(), vec.end());
  ingroup.clear();
  for (int i = fields.indivOffset; i < fields.numIndivs + fields.indivOffset; ++i) {
    if (!std::binary_search(vec.begin(), vec.end(), i)) ingroup.push_back(i);
  }
  std::sort(ingroup.begin(), ingroup.end());
}

std::ostream &vmatch::operator<<(std::ostream &os, const vmatch::vcfrecord &vrec) {
  os << "POS: " << vrec.pos << std::endl;
  os << "REF: P(" << vrec.ref << ")=" << vrec.refFreq << std::endl;
  os << "ALTS: " << std::endl;
  for (auto &e : vrec.altIndivs) {
    os << "\tP(" << e.first << ")=" << vrec.altAF.at(e.first);
    for (auto i : e.second) {
      os << ", " << i;
    }
    os << std::endl;
  }
  return os;
}

vmatch::vcfstream &vmatch::operator>>(vmatch::vcfstream &vstream, vcfrecord &vrec) {
  if (!vstream.getRecord(vrec)) throw std::out_of_range("No more records left.");
  return vstream;
}

void vmatch::vcfstream::printIngroup(std::ostream &os) {
  os << "#";
  for (auto &e : ingroup) {
    os << e << ",";
  }
  os << std::endl;
}