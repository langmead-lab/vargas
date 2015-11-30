/**
 * Ravi Gaddipati
 * November 23, 2015
 * rgaddip1@jhu.edu
 *
 *contains common functions.
 *
 * utils.h
 */

#ifndef VARGAS_UTILS_H
#define VARGAS_UTILS_H

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include "getopt_pp.h"
#include "../gssw/src/gssw.h"
#include "xcoder.h"
#include <boost/algorithm/string/split.hpp>

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


inline void printEncodedVec(std::ostream &os, const std::vector<uint32_t> &vec, vargas::Xcoder &x) {
  if (vec.size() == 0) return;
  os << ':';
  os << x.compressAndEncode(vec);
}


#endif //VARGAS_UTILS_H
