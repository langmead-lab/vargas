//
// Created by gaddra on 9/12/15.
//

#include <iostream>
#include <vector>
#include <sstream>
#include "utils.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  /** Split string with delim, return a vector **/
  std::stringstream ss(s);
  std::string item;
  elems = *new std::vector<std::string>(0);

  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  if (s.at(s.length() - 1) == ',') elems.push_back("");
  return elems;
}
