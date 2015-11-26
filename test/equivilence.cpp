//
// Created by gaddra on 11/25/15.
//

#include "equivilence.h"

bool fileEquiv(std::string afile, std::string bfile, bool orderImp) {
  std::ifstream a(afile);
  std::ifstream b(bfile);
  if (!a.good() || !b.good()) throw std::invalid_argument("Invalid file(s).");
  return fileEquiv(a, b, orderImp);
}

bool fileEquiv(std::string afile, std::istream &b, bool orderImp) {
  std::ifstream s(afile);
  if (!s.good()) throw std::invalid_argument("Invalid file(s).");
  return fileEquiv(s, b, orderImp);
}

bool fileEquiv(std::istream &a, std::string bfile, bool orderImp) {
  std::ifstream s(bfile);
  if (!s.good()) throw std::invalid_argument("Invalid file(s).");
  return fileEquiv(a, s, orderImp);
}

bool fileEquiv(std::istream &a, std::istream &b, bool orderImp) {
  std::string line;
  std::vector<std::string> vec;
  std::set<std::string> set;

  while (std::getline(a, line)) {
    if (line.substr(0, 2) != "##") {
      if (orderImp) vec.push_back(line);
      else set.insert(line);
    }
  }
  int counter = 0;
  while (std::getline(b, line)) {
    if (line.substr(0, 2) != "##") {
      if (orderImp) {
        if (line != vec[counter]) {
          std::cerr << "Line " << counter << " mismatch." << std::endl;
          std::cerr << "\tA: " << vec[counter] << std::endl;
          std::cerr << "\tB: " << line << std::endl;
          return false;
        }
      }
      else {
        if (set.erase(line) != 1) {
          std::cerr << "Line not found in A: " << line << std::endl;
          return false;
        }
      }
      counter++;
    }
  }

  if (orderImp && vec.size() != counter) {
    std::cerr << "Line count mismatch. A: " << vec.size() << ", B: " << counter << std::endl;
    return false;
  }
  else if (!orderImp && set.size() != 0) {
    std::cerr << "File A has " << set.size() << " extra lines." << std::endl;
    return false;
  }

  return true;
}
