//
// Created by gaddra on 6/2/16.
//

#include <vector>
#include <bitset>
#include <sstream>
#include "doctest/doctest.h"

#ifndef VARGAS_DYN_BITSET_H
#define VARGAS_DYN_BITSET_H

/**
 * Groups together multiple fixed bitsets to emulate a dynamic bitset.
 * A vector of std::bitset<core_size> maintains the information.
 * The data is a set with a vector.
 */
template<unsigned int core_size>
class dyn_bitset {

 public:
  dyn_bitset() { }

  dyn_bitset(const dyn_bitset &db) {
    _bitset = db._bitset;
    _right_pad = db._right_pad;
  }

  /**
   * Create a dyn_bitset from a vector. Each bit is
   * set according to the truth of each vector element.
   * @param vector<T>
   */
  template<typename T>
  dyn_bitset(const std::vector<T> &vec) {
    _init_from_vec(vec);
  }

  /**
   * Create a dyn_bitset from a vector. Each bit is
   * set according to the truth of each vector element.
   * @param vector<T>
   */
  template<typename T>
  void set(const std::vector<T> &vec) {
    _init_from_vec(vec);
  }

  /**
   * @return right padding of the rightmost element in the vector
   */
  int right_pad() const { return _right_pad; }

  /**
   * @return number of bits used.
   */
  size_t size() const { return (_bitset.size() * core_size) - _right_pad; }

  /**
   * Set a single bit, default true
   * @param bit bit index
   * @param val value to set bit to.
   */
  void set(const int bit,
           bool val = true) {
    if (bit > size() || bit < 0) throw std::range_error("Index out of bounds.");
    _bitset[bit / 32][bit % 32] = val;
  }

  /**
   * Flips the specified bit.
   * @param bit bit index
   */
  void flip(const int bit) {
    if (bit > size() || bit < 0) throw std::range_error("Index out of bounds.");
    _bitset[bit / 32][bit % 32] = !_bitset[bit / 32][bit % 32];
  }

  /**
   * Value of a single bit
   * @param bit bit index
   */
  bool at(const int bit) const {
    if (bit > size() || bit < 0) throw std::range_error("Index out of bounds.");
    return _bitset[bit / 32][bit % 32];
  }

  /**
   * Const ref to raw data container
   * @return vector of bitsets.
   */
  const std::vector<std::bitset<core_size>> &bitset() const {
    return _bitset;
  }

  bool operator==(const dyn_bitset &db) const {
    if (_bitset.size() != db._bitset.size() || _right_pad != db._right_pad) return false;
    for (size_t i = 0; i < _bitset.size(); ++i) {
      if (_bitset[i] != db._bitset[i]) return false;
    }
    return true;
  }

  /**
   * operator returns true if there is a common bit set.
   * Returns false if there is a dimension mismatch.
   * e.g:
   * 0101 && 1010 = false
   * 0010 && 1011 = true
   */
  bool operator&&(const dyn_bitset &db) const {
    if (_bitset.size() != db._bitset.size() || _right_pad != db._right_pad) return false;
    for (size_t i = 0; i < _bitset.size(); ++i) {
      std::bitset<core_size> t = _bitset[i] & db._bitset[i];
      if (t.any()) return true;
    }
    return false;
  }

  std::string print() const {
    std::stringstream ss;
    for (size_t c = 0; c < _bitset.size(); ++c) {
      for (size_t i = 0; i < core_size; ++i) {
        if (c != _bitset.size() - 1 || i < core_size - _right_pad) {
          ss << (_bitset[c][i] == true ? "1" : "0");
        }
      }
    }
    return ss.str();
  }

 protected:

  template<typename T>
  void _init_from_vec(const std::vector<T> &vec) {
    _bitset.clear();
    _bitset.resize((vec.size() / 32) + 1);
    for (size_t i = 0; i < vec.size(); ++i) {
      const size_t &set_num = i / 32;
      const size_t &bit_offset = i % 32;
      if (vec[i]) _bitset[set_num][bit_offset] = true;
    }
    _right_pad = (_bitset.size() * core_size) - vec.size();
  }

 private:
  std::vector<std::bitset<core_size>> _bitset;
  int _right_pad = 0;
};

TEST_CASE ("Dynamic Bitset") {
  std::vector<bool> a_bool = {0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1,
                              0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
                              1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0,
                              0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0};
  std::vector<bool> b1_bool = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::vector<bool> b2_bool = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0,
                               0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0};
  dyn_bitset<32> a(a_bool),
      b1(std::vector<bool>(62, false)),
      b2(b2_bool);

      CHECK(a.size() == 62);
      CHECK(a.right_pad() == 2);
      CHECK((a && b1) == false);
      CHECK((a && b2) == true);

      CHECK(b1.at(3) == 0);
  b1.set(3);
      CHECK(b1.at(3) == 1);
      CHECK((a && b1) == true);

      CHECK(a.at(0) == 0);
      CHECK(a.at(61) == 0);
      CHECK(a.at(57) == 1);
      CHECK_THROWS(a.at(-1) == 0);
      CHECK_THROWS (a.at(100) == 0);

  std::stringstream ss;
  std::vector<bool> p_bool = {0, 1, 0, 1, 1, 1, 1};
  dyn_bitset<8> p(p_bool);
      CHECK(p.print() == "0101111");
}

#endif //VARGAS_DYN_BITSET_H
