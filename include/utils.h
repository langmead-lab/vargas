/**
 * Ravi Gaddipati
 * November 23, 2015
 * rgaddip1@jhu.edu
 *
 *contains common functions.
 *
 * utils.h
 */

#ifndef VMATCH_UTILS_H
#define VMATCH_UTILS_H

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include "getopt_pp.h"
#include "../gssw/src/gssw.h"

/// <summary>
/// Splits the specified string, resets elems and returns with split string.
/// </summary>
/// <param name="s">The string</param>
/// <param name="delim">The delimiter</param>
/// <param name="elems">Vector to store results in. Vector is replaced!</param>
/// <returns>Vector of split string.</returns>
std::vector<std::string> split(const std::string &s, char delim);

/// <summary>
/// Splits the specified string, resets vec and returns with split string.
/// </summary>
/// <param name="s">The string</param>
/// <param name="delim">The delimiter</param>
/// <param name="vec">Vector to store results in. Vector is cleared!</param>
/// <returns>Vector of split string.</returns>
void split(const std::string &s, char delim, std::vector<std::string> &vec);

template<typename T>
void printVec(std::ostream &os, const std::vector<T> &vec) {
  os << '[';
  for (int i = 0; i < vec.size(); ++i) {
    os << vec[i];
    if (i < vec.size() - 1) os << ',';
  }
  os << ']' << std::endl;
}


#endif //VMATCH_UTILS_H
