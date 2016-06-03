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

#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include "xcoder.h"
#include "doctest/doctest.h"

typedef unsigned char uchar;

/**
 * Converts a character to a numeral representation.
 * @param c character
 * @return numeral representation
 */
__attribute__((always_inline))
inline uchar base_to_num(char c) {
  switch (c) {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      return 4;
  }
}

/**
 * Convert a numeric form to a char, upper case.
 * All ambiguous bases are represented as 'N'
 * @param num numeric form
 * @return char in [A,C,G,T,N]
 */
__attribute__((always_inline))
inline char num_to_base(uchar num) {
  switch (num) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
    default:
      return 'N';
  }
}

/**
 * Convert a character sequence in a numeric sequence
 * @param seq Sequence string
 * @return vector of numerals
 */
__attribute__((always_inline))
inline std::vector<uchar> seq_to_num(const std::string &seq) {
  std::vector<uchar> num(seq.length());
  std::transform(seq.begin(), seq.end(), num.begin(), base_to_num);
  return num;
}
TEST_CASE ("Sequence to Numeric") {
  std::vector<uchar> a = seq_to_num("ACGTN");
      REQUIRE(a.size() == 5);
      CHECK(a[0] == 0);
      CHECK(a[1] == 1);
      CHECK(a[2] == 2);
      CHECK(a[3] == 3);
      CHECK(a[4] == 4);
}

/**
 * Convert a numeric vector to a sequence of bases.
 * @param num Numeric vector
 * @return sequence string, Sigma={A,G,T,C,N}
 */
__attribute__((always_inline))
inline std::string num_to_seq(const std::vector<uchar> &num) {
  std::stringstream builder;
  for (auto &n : num) {
    builder << num_to_base(n);
  }
  return builder.str();
}
TEST_CASE ("Numeric to Sequence") {
  std::string a = num_to_seq({0, 1, 2, 3, 4});
      REQUIRE(a.length() == 5);
      CHECK(a[0] == 'A');
      CHECK(a[1] == 'C');
      CHECK(a[2] == 'G');
      CHECK(a[3] == 'T');
      CHECK(a[4] == 'N');
}


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

std::string getLastLine(std::ifstream& in);


#endif //VARGAS_UTILS_H
