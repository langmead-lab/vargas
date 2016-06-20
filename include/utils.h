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

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include "doctest/doctest.h"

typedef unsigned char uchar;
enum Base: uchar { A = 0, C = 1, G = 2, T = 3, N = 4 };

/**
 * Converts a character to a numeral representation.
 * @param c character
 * @return numeral representation
 */
__attribute__((always_inline))
inline Base base_to_num(char c) {
  switch (c) {
    case 'A':
    case 'a':
      return Base::A;
    case 'C':
    case 'c':
      return Base::C;
    case 'G':
    case 'g':
      return Base::G;
    case 'T':
    case 't':
      return Base::T;
    default:
      return Base::N;
  }
}

/**
 * Convert a numeric form to a char, upper case.
 * All ambiguous bases are represented as 'N'
 * @param num numeric form
 * @return char in [A,C,G,T,N]
 */
__attribute__((always_inline))
inline char num_to_base(Base num) {
  switch (num) {
    case Base::A:
      return 'A';
    case Base::C:
      return 'C';
    case Base::G:
      return 'G';
    case Base::T:
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
inline std::vector<Base> seq_to_num(const std::string &seq) {
  std::vector<Base> num(seq.length());
  std::transform(seq.begin(), seq.end(), num.begin(), base_to_num);
  return num;
}

/**
 * Convert a numeric vector to a sequence of bases.
 * @param num Numeric vector
 * @return sequence string, Sigma={A,G,T,C,N}
 */
__attribute__((always_inline))
inline std::string num_to_seq(const std::vector<Base> &num) {
  std::stringstream builder;
  for (auto &n : num) {
    builder << num_to_base(n);
  }
  return builder.str();
}

TEST_CASE ("Sequence to Numeric") {
  std::vector<Base> a = seq_to_num("ACGTN");
      REQUIRE(a.size() == 5);
      CHECK(a[0] == Base::A);
      CHECK(a[1] == Base::C);
      CHECK(a[2] == Base::G);
      CHECK(a[3] == Base::T);
      CHECK(a[4] == Base::N);
}

TEST_CASE ("Numeric to Sequence") {
  std::string a = num_to_seq({Base::A, Base::C, Base::G, Base::T, Base::N});
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


inline bool file_exists(std::string filename) {
  std::ifstream f(filename);
  return f.good();
}

inline char rand_base() {
  switch (rand() % 5) {
    case 0:
      return 'A';
    case 1:
      return 'T';
    case 2:
      return 'C';
    case 3:
      return 'G';
    default:
      return 'N';
  }
}

/**
 * Edit distance between strings. Taken from:
 * https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C.2B.2B
 */
int levenshtein_distance(const std::string &s1, const std::string &s2);


#endif //VARGAS_UTILS_H
