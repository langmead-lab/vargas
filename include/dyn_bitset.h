/**
 * @file
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date June 05, 2016
 *
 * @brief
 * Emulates a dynamic bitset by utilizing a vector of fixed size bitsets.
 */

#include <vector>
#include <bitset>
#include <sstream>
#include "doctest.h"

#ifndef VARGAS_DYN_BITSET_H
#define VARGAS_DYN_BITSET_H

/**
 * @brief Dynamic bitset backed by fixed bitsets.
 * @details
 * Groups together multiple fixed bitsets to emulate a dynamic bitset.
 * A vector of std::bitset<core_size> maintains the information.
 */
template<unsigned int core_size>
class dyn_bitset {

  public:
    dyn_bitset() = default;
    /**
     * @brief
     * Initilize a bitset of length len, all set to val.
     * @param len bitset length
     * @param val true/false
     */
    dyn_bitset(size_t len,
               bool val = false) : _bitset(std::vector<std::bitset<core_size>>((len / core_size) + 1)) {
        for (std::bitset<core_size> &n : _bitset) {
            if (val) n.set();
            else n.reset();
        }
        _right_pad = (_bitset.size() * core_size) - len;
    };

    /**
     * @brief
     * Create a dyn_bitset from a vector. Each bit is
     * set according to the truth of each vector element.
     * @param vec vector<T>
     */
    template<typename T>
    dyn_bitset(const std::vector<T> &vec) {
        _init_from_vec(vec);
    }

    /**
     * @brief
     * Create a dyn_bitset from a vector. Each bit is
     * set according to the truth of each vector element.
     * @param vec vector<T>
     */
    template<typename T>
    void set(const std::vector<T> &vec) {
        _init_from_vec(vec);
    }

    /**
     * @return right padding of the rightmost element in the vector
     */
    size_t right_pad() const { return _right_pad; }

    /**
     * @return number of bits used.
     */
    inline size_t size() const { return (_bitset.size() * core_size) - _right_pad; }

    /**
     * @brief
     * Set all bits to true.
     */
    void set() {
        for (std::bitset<core_size> &bs : _bitset) bs.set();
    }

    /**
     * @brief
     * Set all bits to false.
     */
    void reset() {
        for (std::bitset<core_size> &bs : _bitset) bs.reset();
    }

    /**
     * @brief
     * Set a single bit, default true
     * @param bit bit index
     * @param val value to set bit to.
     * @throws std::range_error bit is out of range of bitset size
     */
    void set(const size_t bit,
             bool val = true) {
        if (bit > size() || bit < 0) throw std::range_error("Index out of bounds.");
        _bitset[bit / core_size][bit % core_size] = val;
    }

    /**
     * @return true if any bits are set
     */
    bool any() const {
        for (auto &b : _bitset) {
            if (b.any()) return true;
        }
        return false;
    }

    /**
     * @brief
     * Flips the specified bit.
     * @param bit bit index
     * @throws std::range_error bit is out of range
     */
    void flip(const size_t bit) {
        if (bit > size() || bit < 0) throw std::range_error("Index out of bounds.");
        _bitset[bit / core_size][bit % core_size] = !_bitset[bit / core_size][bit % core_size];
    }

    /**
     * @brief
     * Value of a single bit
     * @param bit bit index
     * @throws std::range_error bit is out of range
     */
    bool at(const size_t bit) const {
        if (bit > size() || bit < 0) throw std::range_error("Index out of bounds.");
        return _bitset[bit / core_size][bit % core_size];
    }

    /**
     * @brief
     * Alias for at(). Does not allow setting.
     * @param bit index of bit to get
     */
    bool operator[](const size_t bit) const {
        return at(bit);
    }

    /**
     * @brief
     * Flips every bit (negation)
     * @return negated bitset
     */
    dyn_bitset<core_size> operator~() const {
        dyn_bitset<core_size> ret = *this;
        for (size_t i = 0; i < size(); ++i) {
            ret.flip(i);
        }
        return ret;
    }

    /**
     * @brief
     * Add a bit.
     * @param val 0/1 bit
     */
    void push_back(const bool val) {
        if (_right_pad == 0) {
            _bitset.push_back(std::bitset<core_size>());
            _right_pad = core_size;
        }
        _bitset[size() / core_size][size() % core_size] = val;
        --_right_pad;
    }

    /**
     * brief
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

    bool operator!=(const dyn_bitset &db) const { return !operator==(db); }

    /**
     * @brief
     * operator returns true if there is a common bit set.
     * Returns false if there is a dimension mismatch.
     * @details
     * e.g: \n
     * 0101 && 1010 = false \n
     * 0010 && 1011 = true
     * @param db other bitset
     * @throws std::invalid_argument bitsets are different sizes
     */
    bool operator&&(const dyn_bitset &db) const {
        if (size() != db.size())
            throw std::invalid_argument("Incompatible dimension");
        for (size_t i = 0; i < _bitset.size(); ++i) {
            std::bitset<core_size> t = _bitset[i] & db._bitset[i];
            if (t.any()) {
                return true;
            }
        }
        return false;
    }

    /**
     * Bitwise AND
     * @param other
     * @return b1 & b2
     * @throws std::range_error Bitsets are incompatible dimensions
     */
    dyn_bitset<core_size> operator&(const dyn_bitset<core_size> &other) {
        //TODO done slowly
        //TODO add tests
        if (size() != other.size())
            throw std::range_error("Incompatible dimension");
        dyn_bitset<core_size> ret(size());
        for (size_t i = 0; i < size(); ++i) {
            if (at(i) && other.at(i)) {
                ret.set(i);
            }
        }
        return ret;
    }

    /**
     * Bitwise OR
     * @param other
     * @return b1 | b2
     * @throws std::range_error Bitsets are incompatible dimensions
 */
    dyn_bitset<core_size> operator|(const dyn_bitset<core_size> &other) {
        //TODO done slowly
        //TODO add tests
        if (size() != other.size())
            throw std::range_error("Incompatible dimension");
        dyn_bitset<core_size> ret(size());
        for (size_t i = 0; i < size(); ++i) {
            if (at(i) || other.at(i)) {
                ret.set(i);
            }
        }
        return ret;
    }

    /**
     * @param pop init via assignment from vector
     */
    template<typename T>
    dyn_bitset &operator=(const std::vector<T> &pop) {
        _init_from_vec<T>(pop);
        return *this;
    }

    /**
     * @brief
     * Count the number of bits set.
     * @return number of bits set.
     */
    size_t count() const {
        size_t count = 0;
        for (size_t i = 0; i < size(); ++i) {
            if (at(i)) ++count;
        }
        return count;
    }


    /**
     * @return string of 0's and 1's
     */
    std::string to_string() const {
        std::ostringstream ss;
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

    /**
     * @brief
     * Sets a bit if the corresponding vector element tests true.
     * @param vec vector of values to test
     */
    template<typename T>
    void _init_from_vec(const std::vector<T> &vec) {
        _bitset.clear();
        _bitset.resize((vec.size() / core_size) + 1);
        for (size_t i = 0; i < vec.size(); ++i) {
            const size_t &set_num = i / core_size;
            const size_t &bit_offset = i % core_size;
            if (vec[i]) _bitset[set_num][bit_offset] = true;
        }
        _right_pad = (_bitset.size() * core_size) - vec.size();
    }

  private:
    std::vector<std::bitset<core_size>> _bitset;
    size_t _right_pad = 0;
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

    a.push_back(0);
    CHECK(a.size() == 63);
    CHECK(a.at(62) == 0);
    a.push_back(1);
    CHECK(a.size() == 64);
    CHECK(a.at(63) == 1);
    a.push_back(0);
    CHECK(a.size() == 65);
    CHECK(a.at(64) == 0);

    std::ostringstream ss;
    std::vector<bool> p_bool = {0, 1, 0, 1, 1, 1, 1};
    dyn_bitset<8> p(p_bool);
    CHECK(p.to_string() == "0101111");
}

#endif //VARGAS_DYN_BITSET_H
