/**
 * Ravi Gaddipati
 * November 23, 2015
 * rgaddip1@jhu.edu
 *
 *contains common functions.
 *
 * utils.cpp
 */


#include "../include/utils.h"


std::vector<std::string> split(const std::string &s, char delim) {
  /** Split string with delim, return a vector **/
  std::vector<std::string> newElems(0);
  split(s, delim, newElems);
  return newElems;
}


inline void split(const std::string &s, char delim, std::vector<std::string> &vec) {
  /** Split string with delim, return a vector **/
  std::stringstream ss(s);
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