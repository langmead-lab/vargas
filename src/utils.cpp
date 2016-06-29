/**
 * @author Ravi Gaddipati
 * @date November 23, 2015
 * rgaddip1@jhu.edu
 *
 * @brief
 * Contains common functions.
 *
 * @file
 */


#include <sstream>
#include "utils.h"

std::vector<std::string> split(const std::string &s, char delim) {
  /** Split string with delim, return a vector **/
  std::vector<std::string> newElems(0);
  split(s, delim, newElems);
  return newElems;
}


inline void split(const std::string &s, char delim, std::vector<std::string> &vec) {
  /** Split string with delim, return a vector **/
  std::istringstream ss(s);
  std::string item;
  vec.clear();


  if (s.length() == 0) {
    return;
  }
  else if (s.length() == 1 && s.at(0) != delim) {
    vec.push_back(s.substr(0, 1));
    return;
  }

  while (std::getline(ss, item, delim)) {
    if (item.size() != 0) vec.push_back(item);
  }

}

std::string getLastLine(std::ifstream& in) {
  std::string line, ret = "";
  while (in >> std::ws && std::getline(in, line)) {
    if (line.at(0) != '#') ret = line;
  }
  return ret;
}

int levenshtein_distance(const std::string &s1, const std::string &s2) {
  // To change the type this function manipulates and returns, change
  // the return type and the types of the two variables below.
  int s1len = s1.size();
  int s2len = s2.size();

  auto column_start = (decltype(s1len)) 1;

  auto column = new decltype(s1len)[s1len + 1];
  std::iota(column + column_start, column + s1len + 1, column_start);

  for (auto x = column_start; x <= s2len; x++) {
    column[0] = x;
    auto last_diagonal = x - column_start;
    for (auto y = column_start; y <= s1len; y++) {
      auto old_diagonal = column[y];
      auto possibilities = {
          column[y] + 1,
          column[y - 1] + 1,
          last_diagonal + (s1[y - 1] == s2[x - 1] ? 0 : 1)
      };
      column[y] = std::min(possibilities);
      last_diagonal = old_diagonal;
    }
  }
  auto result = column[s1len];
  delete[] column;
  return result;
}