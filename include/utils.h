/**
 * Ravi Gaddipati
 * June 26, 2016
 * rgaddip1@jhu.edu
 *
 * @brief
 * Contains common functions.
 *
 * @copyright
 * Distributed under the MIT Software License.
 * See accompanying LICENSE or https://opensource.org/licenses/MIT
 *
 * @file
 */

#ifndef VARGAS_UTILS_H
#define VARGAS_UTILS_H

#define RG_UTIL_INCLUDE_DOCTESET 1

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
#if defined(__GNUC__)
#define __RG_STRONG_INLINE__ __attribute__((always_inline)) inline
#else
#define __RG_STRONG_INLINE__ inline
#endif
#else
#pragma message("No strong inlining")
#define __RG_STRONG_INLINE__ inline
#endif

#include <array>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <random>
#include <type_traits>
#include <cassert>
#include <memory>
#include <iostream>
#include <stdexcept>

namespace rg {

  using pos_t = uint32_t;

/**
 * @enum Base
 * @brief
 * Maps Base characters to integers.
 * @warning Numeric assignments should not be changed
 */
  enum Base: unsigned char { N = 0, A = 1, C = 2, G = 3, T = 4 };

/**
 * @brief
 * Converts a character to a numeral representation.
 * @param c character
 * @return numeral representation
 */
  __RG_STRONG_INLINE__
  Base base_to_num(char c) {
/* Generated with:
lut = ["Base::N"] * 128;
lut[ord('A')] = lut[ord('a')] = "Base::A"
lut[ord('C')] = lut[ord('c')] = "Base::C"
lut[ord('G')] = lut[ord('g')] = "Base::G"
lut[ord('T')] = lut[ord('t')] = "Base::T"
print("        " + ", ".join(lut))
*/
    static constexpr std::array<Base, 128> lut {
        Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::A, Base::N, Base::C, Base::N, Base::N, Base::N, Base::G, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::T, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::A, Base::N, Base::C, Base::N, Base::N, Base::N, Base::G, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::T, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N, Base::N
    };
    return lut[static_cast<uint8_t>(c)];
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
      // Leaving alone, considering this encoding might change.
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
      std::string ret; ret.reserve(num.size());
      for(const auto c: num) ret += num_to_base(c);
      return ret;
  }

  // TODO: replace with LUT

static constexpr std::array<char, 128> nuc_complement_lut {
    78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 84, 78, 71, 78, 78, 78, 67, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 65, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78
};
  __RG_STRONG_INLINE__
  char complement(const char b) {
// Auto-generated from the following Python code:
/*
lut = ['N'] * 128
lut[ord('A')] = lut[ord('a')] = 'T'
lut[ord('C')] = lut[ord('c')] = 'G'
lut[ord('G')] = lut[ord('g')] = 'C'
lut[ord('T')] = lut[ord('t')] = 'A'
print(", ".join(map(str, map(ord, lut))))
*/
    return nuc_complement_lut[b];
  }

  __RG_STRONG_INLINE__
  Base complement_b(const Base b) {
      switch (b) {
          case Base::A:
              return Base::T;
          case Base::C:
              return Base::G;
          case Base::G:
              return Base::C;
          case Base::T:
              return Base::A;
          default:
              return Base::N;
      }
  }

  //TODO resolve this
  __RG_STRONG_INLINE__
  void reverse_complement_inplace(std::string &seq){
      std::transform(seq.begin(), seq.end(), seq.begin(), complement);
      std::reverse(seq.begin(), seq.end());
      /* // Reverted because executing this was causing a segfault later on
      size_t i;
      for(i = 0; i < seq.size() << 1; ++i) {
          char tmp = seq[seq.size() - i - 1];
          seq[seq.size() - i - 1] = nuc_complement_lut[seq[i]];
          seq[i] = nuc_complement_lut[tmp];
      }
      if(i & 1) seq[i] = nuc_complement_lut[seq[i]];
       */
      // Avoids two passes
  }

  __RG_STRONG_INLINE__
  std::string reverse_complement(const std::string &seq) {
      std::string ret = seq;
      reverse_complement_inplace(ret);
      return ret;
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
      vec = split(s, std::string(1, delim));
  }

/**
 * @brief
 * Splits a string into a vector and guesses the delimiter.
 * @param s string to split
 * @param vec vector to store results in
 */
  inline void split(const std::string &s, std::vector<std::string> &vec) {
      vec = split(s, std::string(1, rg::guess_delim(s)));
  }


/**
 * @brief
 * Opens a file and checks if its valid.
 * @param filename File to check if valid.
 */
  inline bool file_exists(const std::string& filename) {
      std::ifstream f(filename);
      return f.good();
  }

/**
 * @return random base character in [ACGTN].
 */
  __RG_STRONG_INLINE__
  char rand_base() {
      // Avoid use of std::rand
      static std::minstd_rand rng(1337);
      static uint32_t cval = rng();
      static int nleft = 10;
      if(nleft == 0) {
          cval = rng();
          nleft = sizeof(typename std::minstd_rand::result_type) / 3;
      }
      const auto sval = cval % 5;
      cval >>= 3;
      --nleft;
      return "ACGTN"[sval];
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
  typename std::enable_if<std::is_arithmetic<T>::value, std::string>::type
  to_string(T val) {
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

  template<typename T>
  std::string vec_to_str(const std::vector<T> &v, const std::string& sep = ", ") {
      std::stringstream os;
      if (!v.size()) return "";
      os << v[0];
      for (unsigned i = 1; i < v.size(); ++i) {
          os << sep << v[i];
      }
      return os.str();
  }

  template<typename T>
  void print_vec(const std::vector<T> &v, std::ostream &os = std::cerr) {
      os << vec_to_str(v);
  }

  template<typename Object, typename R, typename ...Args>
  struct smart_fun {
      Object & obj;
      R (Object::*fun)(Args...);
      R operator()(Args... args) {
          return (obj.*fun)(args...);
      }
  };


  /**
   * @brief
   * Return a function proxy for class members.
   * @details
   * Usage:
   * Foo foo;
   * auto fn = smart_bind(foo, Foo::bar);
   * Allows stateful function calls.
   */
  template<typename C, typename R, typename ...Args>
  auto smart_bind(C & c, R (C::*fun)(Args...)) -> smart_fun<C, R, Args...> {
      return smart_fun<C, R, Args...>{c, fun};
  }

  /**
   * make_unique, from:
   * https://isocpp.org/files/papers/N3656.txt
   */
  template<class T>
  struct _Unique_if {
      typedef std::unique_ptr<T> _Single_object;
  };

  template<class T>
  struct _Unique_if<T[]> {
      typedef std::unique_ptr<T[]> _Unknown_bound;
  };

  template<class T, size_t N>
  struct _Unique_if<T[N]> {
      typedef void _Known_bound;
  };

  template<class T, class... Args>
  typename _Unique_if<T>::_Single_object
  make_unique(Args &&... args) {
      return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  template<class T>
  typename _Unique_if<T>::_Unknown_bound
  make_unique(size_t n) {
      typedef typename std::remove_extent<T>::type U;
      return std::unique_ptr<T>(new U[n]());
  }

  template<class T, class... Args>
  typename _Unique_if<T>::_Known_bound
  make_unique(Args &&...) = delete;

  struct Deleter {
      void operator()(const void *p) const {
          ::std::free(const_cast<void *>(p));
      }
  };

}

template<typename T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    rg::print_vec(v, os);
    return os;
}

#endif //VARGAS_UTILS_H
