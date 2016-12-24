/**
 * Ravi Gaddipati
 * June 26, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Contains common functions.
 *
 * @file
 */

#ifndef VARGAS_UTILS_H
#define VARGAS_UTILS_H

#if defined __INTEL_COMPILER
#define __RG_UNROLL__
#define RESTRICT restrict
#elif defined __GNUC__
#define __RG_UNROLL__ __attribute__((optimize("unroll-loops")))
#define RESTRICT __restrict__
#else
#define __RG_LIKELY__(x) (x)
#define __RG_UNLIKELY__(x) (x)
#define __RG_UNROLL__
#define RESTRICT
#endif

#ifndef RG_DISABLE_INLINE
#if defined __GNUC__
#define __RG_STRONG_INLINE__ __attribute__((always_inline)) inline
#else
#define __RG_STRONG_INLINE__ inline
#endif
#else
#define __RG_STRONG_INLINE__ inline
#endif

#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <type_traits>
#include <assert.h>
#include "simdpp/simd.h"

namespace rg {

/**
 * @enum Base
 * @brief
 * Maps Base characters to integers.
 */
  enum Base: unsigned char { A = 0, C = 1, G = 2, T = 3, N = 4 };

/**
 * @brief
 * Converts a character to a numeral representation.
 * @param c character
 * @return numeral representation
 */
  __RG_STRONG_INLINE__
  Base base_to_num(char c) {
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
 * @brief
 * Convert a numeric form to a char, upper case.
 * All ambiguous bases are represented as 'N'.
 * @param num numeric form
 * @return char in [A,C,G,T,N]
 */
  __RG_STRONG_INLINE__
  char num_to_base(Base num) {
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
 * @brief
 * Convert a character sequence in a numeric sequence.
 * @param seq Sequence string
 * @return vector of numerals
 */
  __RG_STRONG_INLINE__
  std::vector<Base> seq_to_num(const std::string &seq) {
      std::vector<Base> num(seq.length());
      std::transform(seq.begin(), seq.end(), num.begin(), base_to_num);
      return num;
  }

/**
 * @brief
 * Convert a numeric vector to a sequence of bases.
 * @param num Numeric vector
 * @return sequence string, Sigma={A,G,T,C,N}
 */
  __RG_STRONG_INLINE__
  std::string num_to_seq(const std::vector<Base> &num) {
      std::ostringstream builder;
      for (auto &n : num) {
          builder << num_to_base(n);
      }
      return builder.str();
  }


  /**
   * @brief
   * Apply binary_op to each token in [in_begin,in_end], split by any element in [d_begin, d_end].
   * Adapted from http://tristanbrindle.com/posts/a-quicker-study-on-tokenising/
   * @tparam InIter Type of input iterator
   * @tparam DelimIter Type of split iterator
   * @tparam BinOp Binary Op to apply fo each token
   * @param in_begin Input begin iterator
   * @param in_end Input end iterator
   * @param d_begin Delimiter begin iterator
   * @param d_end Delimiter end iterator
   * @param binary_op operation to apply to each token
   */
  template<typename InIter, typename DelimIter, class BinOp>
  void for_each_token(InIter in_begin, InIter in_end, DelimIter d_begin, DelimIter d_end, BinOp binary_op) {
      while (in_begin != in_end) {
          const auto pos = std::find_first_of(in_begin, in_end, d_begin, d_end);
          binary_op(in_begin, pos);
          if (pos == in_end) break;
          in_begin = std::next(pos);
      }
  }

  /**
   * @brief
   * Split a string into tokens with any delimiter in delims.
   * @param str String to split
   * @param delims List of delimiters to split by
   * @param skip_empty Skip empty elements
   * @return vector of tokens
   */
  std::vector<std::string> split(const std::string &str, const std::string &delims = ",", bool skip_empty = true);

/**
 * @brief
 * Splits a string into a vector given some character delimiter.
 * @param s string to split
 * @param delim split string at delim, discarding the delim
 * @return vector to store results in
 */
  inline std::vector<std::string> split(const std::string &s, char delim) {
      return split(s, std::string(1, delim));
  }

  /**
    * @brief
    * Guesses the delimiter of a line. Checks the following, with the first that succesfully splits
    * the string returned.
    * "<newline><tab>:;,=|/"
    * @param line
    * @return delimiter of the line.
    */
  char guess_delim(const std::string &line);

/**
 * @brief
 * Splits a string into a vector given some character delimiter.
 * @param s string to split
 * @param delim split string at delim, discarding the delim
 * @param vec vector to store results in
 */
  inline void split(const std::string &s, char delim, std::vector<std::string> &vec) {
      vec = std::move(split(s, std::string(1, delim)));
  }

/**
 * @brief
 * Splits a string into a vector and guesses the delimiter.
 * @param s string to split
 * @param vec vector to store results in
 */
  inline void split(const std::string &s, std::vector<std::string> &vec) {
      vec = std::move(split(s, std::string(1, rg::guess_delim(s))));
  }


/**
 * @brief
 * Opens a file and checks if its valid.
 * @param filename File to check if valid.
 */
  inline bool file_exists(std::string filename) {
      std::ifstream f(filename);
      return f.good();
  }

/**
 * @return random base character in [ACGTN].
 */
  __RG_STRONG_INLINE__
  char rand_base() {
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
   * @brief
   * Date in format:
   *    YEAR-MONTH-DAY
   * @return curent date
   */
  std::string current_date();

  /**
   * @tparam T
   * @param start_time
   * @param end_time
   * @return Time between start and end times.
   */
  template<typename T>
  __RG_STRONG_INLINE__
  double chrono_duration(const std::chrono::time_point<T> &start_time, const std::chrono::time_point<T> &end_time) {
      return std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
  }

  /**
   * @tparam T
   * @param start_time
   * @return Time between start time and now.
   */
  template<typename T>
  __RG_STRONG_INLINE__
  double chrono_duration(const std::chrono::time_point<T> &start_time) {
      return chrono_duration(start_time, std::chrono::steady_clock::now());
  }

  /**
   * @brief
   * Checks if a string ends with some substring
   * @param full_string
   * @param ending
   * @return True if full_string ends with ending
   */
  __RG_STRONG_INLINE__
  bool ends_with(std::string const &full_string, std::string const &ending) {
      if (full_string.length() >= ending.length())
          return full_string.compare(full_string.length() - ending.length(), ending.length(), ending) == 0;
      else return false;
  }


  /**
   * @brief
   * Wrapper around std::to_string that allows transparent usage with to_string(std::string)
   * @tparam T Type to convert to string
   * @param val
   * @return
   */
  template<typename T>
  __RG_STRONG_INLINE__
  typename std::enable_if<std::is_arithmetic<T>::value, std::string>::type to_string(T val) {
      return std::to_string(val);
  }

  /**
 * @brief
 * Wrapper around std::to_string that allows transparent usage with to_string(std::string).
 * T = std::string pass through.
 * @tparam T Type to convert to string
 * @param val
 * @return
 */
  __RG_STRONG_INLINE__
  std::string to_string(std::string s) { return s; }

  /**
   * @brief
   * Convert from a string to the respective type.
   * Float types.
   * @tparam T Type to convert to
   * @param s
   * @param ret
   */
  template<typename T>
  __RG_STRONG_INLINE__
  typename std::enable_if<std::is_floating_point<T>::value>::type
  from_string(const std::string &s, T &ret) {
      ret = std::stod(s);
  }

  /**
 * @brief
 * Convert from a string to the respective type.
 * Integer types.
 * @tparam T Type to convert to
 * @param s
 * @param ret
 */
  template<typename T>
  __RG_STRONG_INLINE__
  typename std::enable_if<std::is_integral<T>::value>::type
  from_string(const std::string &s, T &ret) {
      ret = std::stoi(s);
  }

  /**
 * @brief
 * Convert from a string to the respective type.
 * std::string pass through.
 * @tparam T Type to convert to
 * @param s
 * @param ret
 */
  __RG_STRONG_INLINE__
  void from_string(const std::string &s, std::string &ret) { ret = s; }

}

#endif //VARGAS_UTILS_H
