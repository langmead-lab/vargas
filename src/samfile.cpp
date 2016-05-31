/**
 * Ravi Gaddipati
 * April 22, 2016
 * rgaddip1@jhu.edu
 *
 * Interface to work with SAM files.
 *
 * samfile.cpp
 */


#include <vector>
#include "../include/samfile.h"

vargas::Read &vargas::SAMFile::get_read() {
  return currRead;
}
std::string vargas::SAMFile::get_header() const {
  std::stringstream ss;
  for (auto comm : _comments) {
    ss << comm << '\n';
  }
  return ss.str();
}
bool vargas::SAMFile::update_read() {
  std::string line;
  if (!getline(input, line)) return false;
  currAlignment = parseAlignment(line);

  if (currAlignment.flags.reverseComplemented) currRead.read = reverseComplement(currAlignment.seq);
  else currRead.read = currAlignment.seq;

  currRead.end_pos = currAlignment.pos;
  currRead.indiv = -1;
  currRead.indel_err = -1;
  currRead.sub_err = -1;
  currRead.var_bases = -1;
  currRead.var_nodes = -1;
  return true;
}

std::string vargas::SAMFile::reverseComplement(std::string &seq) {
  std::string rc = "";
  for (char &c : seq) {
    if (c == 'A') rc += 'T';
    else if (c == 'T') rc += 'A';
    else if (c == 'G') rc += 'C';
    else if (c == 'C') rc += 'G';
    else if (c == 'a') rc += 't';
    else if (c == 't') rc += 'a';
    else if (c == 'g') rc += 'c';
    else if (c == 'c') rc += 'g';
  }
  return rc;
}

vargas::SAMFile::Alignment vargas::SAMFile::parseAlignment(std::string &line) {
  std::vector<std::string> splitLine = split(line, '\t');
  if (splitLine.size() < 11) throw std::invalid_argument("Invalid alignment line:\n" + line);
  Alignment a;

  a.qname = splitLine[0];
  a.flags = Flag(std::stoul(splitLine[1]));
  a.rname = splitLine[2];
  a.pos = std::stoul(splitLine[3]);
  a.mapq = std::stoi(splitLine[4]);
  a.cigar = splitLine[5];
  a.rnext = std::stoi(splitLine[6]);
  a.pnext = std::stoul(splitLine[7]);
  a.tlen = std::stoi(splitLine[8]);
  a.seq = splitLine[9];
  a.qual = splitLine[10];

  // Load optional fields
  for (int i = 11; i < splitLine.size(); ++i) {
    a.optionalFields[splitLine[i].substr(0, 2)] = splitLine[i].substr(3);
  }

  return a;
}

void vargas::SAMFile::parseHeaderLine(std::string &line) {
  //TODO those if else chains..., can be converted to hardcoded maps
  std::vector<std::string> splitLine;
  std::string tag, val;

  split(line.substr(1), '\t', splitLine);

  if (splitLine[0] == "HD");
  else if (splitLine[0] == "SQ") {
    ReferenceSequence rsq;
    rsq.raw = line;
    for (int i = 1; i < splitLine.size(); ++i) {
      // Parse all the fields
      tag = splitLine[i].substr(0, 2);
      val = splitLine[i].substr(3);
      if (tag == "SN") {
        rsq.name = val;
      } else if (tag == "LN") {
        rsq.length = std::stoi(val);
      } else if (tag == "AS") {
        rsq.genAsm = val;
      } else if (tag == "M5") {
        rsq.md5 = val;
      } else if (tag == "SP") {
        rsq.species = val;
      } else if (tag == "UR") {
        rsq.ur = val;
      }
    }
    refSeqs[rsq.name] = rsq; // each seq name should be unique as per SAM spec
  }

  else if (splitLine[0] == "RG") {
    ReadGroup rg;
    rg.raw = line;
    for (int i = 1; i < splitLine.size(); ++i) {
      // Parse all the fields
      tag = splitLine[i].substr(0, 2);
      val = splitLine[i].substr(3);
      if (tag == "ID") {
        rg.id = val;
      } else if (tag == "CN") {
        rg.name = std::stoi(val);
      } else if (tag == "DS") {
        rg.desc = val;
      } else if (tag == "DT") {
        rg.date = val;
      } else if (tag == "FO") {
        rg.flowOrder = val;
      } else if (tag == "KS") {
        rg.keySeq = val;
      } else if (tag == "LB") {
        rg.lib = val;
      } else if (tag == "PG") {
        rg.prog = val;
      } else if (tag == "PI") {
        rg.median = val;
      } else if (tag == "PL") {
        rg.platform = val;
      } else if (tag == "PM") {
        rg.model = val;
      } else if (tag == "PU") {
        rg.unit = val;
      } else if (tag == "SM") {
        rg.sample = val;
      }
    }
    readGroups[rg.id] = rg; // each id should be unique as per SAM spec
  }
  else if (splitLine[0] == "PG") {
    Program pg;
    pg.raw = line;

    for (int i = 1; i < splitLine.size(); ++i) {
      // Parse all the fields
      tag = splitLine[i].substr(0, 2);
      val = splitLine[i].substr(3);
      if (tag == "ID") {
        pg.id = val;
      } else if (tag == "PN") {
        pg.name = std::stoi(val);
      } else if (tag == "CL") {
        pg.commandLine = val;
      } else if (tag == "PP") {
        pg.prevPrg = val;
      } else if (tag == "DS") {
        pg.desc = val;
      } else if (tag == "VN") {
        pg.version = val;
      }
    }

    programs[pg.id] = pg;
  }
  else {
    _comments.push_back(line);
  }
}

void vargas::SAMFile::init() {
  std::string line;
  while (getline(input, line)) {
    if (line.at(0) == '@') parseHeaderLine(line);
    else {
      currAlignment = parseAlignment(line);
      break;
    }
  }
}