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

#ifdef __GNUC__
// expr x is likely to evaluate true
#define __RG_LIKELY__(x) __builtin_expect((x),1)
// expr x is unlikely to evaluate true
#define __RG_UNLIKELY__(x) __builtin_expect((x),0)
#define __RG_UNROLL__ __attribute__((optimize("unroll-loops")))
#define __RG_STRONG_INLINE__ __attribute__((always_inline)) inline
#else
#define __RG_LIKELY__(x) (x)
#define __RG_UNLIKELY__(x) (x)
#define __RG_UNROLL__
#define __RG_STRONG_INLINE__ inline
#endif

#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <chrono>
#include "simdpp/simd.h"

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
 * Splits a string into a vector given some character delimiter.
 * @param s string to split
 * @param delim split string at delim, discarding the delim
 * @return vector to store results in
 */
std::vector<std::string> split(const std::string &s, char delim);

/**
 * @brief
 * Splits a string into a vector given some character delimiter.
 * @param s string to split
 * @param delim split string at delim, discarding the delim
 * @param vec vector to store results in
 */
void split(const std::string &s,
           char delim, std::vector<std::string> &vec);

/**
 * @brief
 * Splits a string into a vector and guesses the delimiter.
 * @param s string to split
 * @param vec vector to store results in
 */
void split(const std::string &s, std::vector<std::string> &vec);


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
 * Edit distance between strings. Taken from: \n
 * https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C.2B.2B
 * @param s1 sequence a
 * @param s2 sequence b
 */
int levenshtein_distance(const std::string &s1,
                         const std::string &s2);

/**
 * @brief
 * Guesses the delimiter of a line. Checks the following, with the first that succesfully splits
 * the string returned.
 * "<newline><tab>:;,=|/"
 * @param line
 * @return delimiter of the line.
 */
__RG_STRONG_INLINE__
char guess_delim(const std::string &line) {
    static const std::string options = "\n\t:;,=|/";
    for (const char d : options) {
        if (split(line, d).size() > 1) return d;
    }
    throw std::logic_error("Unable to determine delimiter in line: " + line);
}

std::string current_date();

/**
 * @brief
 * Extract the i'th element from a vector. No range checking is done.
 * @param i index of element
 * @param vec vector to extract from
 */
__RG_STRONG_INLINE__
uint8_t extract(uint8_t i,
                const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &vec) {
    return ((uint8_t *) &vec)[i];
}

/**
 * @brief
 * Insert into the i'th element from a vector. No range checking is done.
 * @param elem element to insert
 * @param i index of element
 * @param vec vector to insert in
 */
__RG_STRONG_INLINE__
void insert(uint8_t elem,
            uint8_t i,
            const simdpp::uint8<SIMDPP_FAST_INT8_SIZE> &vec) {
    ((uint8_t *) &vec)[i] = elem;
}


template<typename T>
__RG_STRONG_INLINE__
double chrono_duration(const std::chrono::time_point<T> &start_time) {
    return std::chrono::duration_cast<std::chrono::duration<double>>
        (std::chrono::steady_clock::now() - start_time).count();
}

template<typename T>
__RG_STRONG_INLINE__
double chrono_duration(const std::chrono::time_point<T> &start_time, const std::chrono::time_point<T> &end) {
    return std::chrono::duration_cast<std::chrono::duration<double>>
        (end - start_time).count();
}

__RG_STRONG_INLINE__
bool ends_with(std::string const &fullString,
               std::string const &ending) {
    if (fullString.length() >= ending.length())
        return 0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending);
    else return false;
}

#endif //VARGAS_UTILS_H
