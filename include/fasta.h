/**
 * @file
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date May 26, 2016
 *
 * @brief
 * FASTAFile provides an interface to a FASTA formatted file through htslib.
 * @details
 * An index is created for the opened file if it does not exist.
 */

#ifndef VARGAS_FASTA_H
#define VARGAS_FASTA_H

#include <string>
#include <fstream>
#include <iostream>
#include "utils.h"
#include "doctest.h"
#include "htslib/faidx.h"


namespace Vargas {

  /**
   * @brief
   * FASTA file writer.
   * @details
   * Usage: \n
   * @code{.cpp}
   * #include "fasta.h"
   *
   * {
   *     Vargas::ofasta out("output.fa");
   *     out.char_per_line(5);
   *
   *     out.write("seq1", "ACGTCACCT");
   *     out.write("seq2", "ACGT");
   *
   *     // output.fa:
   *     // >seq1
   *     // ACGTC
   *     // ACCT
   *     // >seq2
   *     // ACGT
   * }
   *
   * {
   *    Vargas ofasta out();
   *    out.char_per_line(3);
   *    out.write("seq1", "ACGTCTCTC");
   *
   *    // stdout:
   *    // >seq1\nACG\nTCT\nCTC\n
   * }
   * @endcode
   */
  class ofasta {
    public:

      /**
       * @brief
       * Output to stdout
       */
      ofasta() : _use_stdio(true) { }

      /**
       * @param file_name output file
       */
      ofasta(std::string file_name) {
          open(file_name);
      }

      /**
       * @brief open a file for output
       * @param file_name
       * @throws std::invalid_argument Error opening a file
       */
      void open(std::string file_name) {
          close();
          if (file_name.length() == 0) {
              _use_stdio = true;
          }
          else {
              _use_stdio = false;
              _o.open(file_name);
              if (!_o.good()) throw std::invalid_argument("Error opening file \"" + file_name + "\"");
          }
      }

      void close() {
          _o.close();
      }

      /**
       * @brief
       * Write a sequence
       * @param name sequence name and any meta information.
       * @param sequence sequence string
       */
      void write(const std::string &name,
                 const std::string &sequence) {
          (_use_stdio ? std::cout : _o) << '>' << name << '\n';
          size_t pos = 0;
          while (pos < sequence.length()) {
              (_use_stdio ? std::cout : _o) << sequence.substr(pos, _char_per_line) << '\n';
              pos += _char_per_line;
          }
      }

      /**
       * @brief
       * Set the number of characters per line.
       * @param len maximum sequence line length
       */
      void char_per_line(int len) {
          if (len > 0) _char_per_line = len;
      }


    private:
      int _char_per_line = 70;
      std::ofstream _o;
      bool _use_stdio;
  };

  /**
   * @brief
   * Provides an interface for a FASTA File. An index is built if one does not
   * already exist.
   * @details
   * Usage:\n
   * @code{.cpp}
   * #include "fasta.h"
   *
   * Vargas::ifasta in("hs37d5.fa"); // input.fa.fai is also loaded, made if it doesn't exist
   * std::vector<std>>string> names = in.sequence_names(); // "1", "2", "3" ... "22", "X", "Y", ...
   * in.num_seq(); // 86 sequences
   *
   * std::string chr22 = in.seq("22"); // All of chromosome 22
   * std::string chr22_0_1000 = in.subseq("22", 0, 1000); // subsequence of a record
   * @endcode
   */
  class ifasta {

    public:

      ifasta() { }

      /**
       * @brief
       * Load a given FASTA and its index, create one if none exists.
       * @param file filename
       */
      ifasta(std::string file) : _file_name(file) {
          open(file);
      }

      ~ifasta() {
          close();
      }

      /**
       * @brief
       * Close any opened file and flush any outputs.
       */
      void close() {
          if (_index) fai_destroy(_index);
          _index = nullptr;
      }

      /**
       * @brief
       * Open a specified FASTA file and make an index.
       * @param file_name filename
       * @return -1 on index build error, -2 on open error, 0 otherwise
       */
      int open(const std::string &file_name) {
          // Check if a Fasta index exists. If it doesn't build it.
          if (!file_exists(file_name + ".fai")) {
              if (fai_build(file_name.c_str()) != 0) {
                  return -1;
              }
          }
          _index = fai_load(file_name.c_str());
          if (!_index) return -2;
          _file_name = file_name;

          _seq_names.clear();
          for (size_t i = 0; i < num_seq(); ++i) {
              _seq_names.push_back(seq_name(i));
          }
          return 0;
      }

      /**
       * @return opened file name.
       */
      std::string file() {
          return _file_name;
      }

      /**
       * @brief
       * Check the number of sequences in the index. In read mode,
       * get the number of buffered sequences.
       * @return number of sequences in the FASTA file
       */
      size_t num_seq() const {
          return faidx_nseq(_index);
      }

      /**
       * @brief
       * Get the sequence with a given name.
       * @param name sequence name
       * @return sequence
       */
      std::string seq(const std::string &name) const {
          return std::string(subseq(name, 0, faidx_seq_len(_index, name.c_str())));
      }

      /**
       * @brief
       * Return a subsequence of a FASTA sequence, 0 based indexing.
       * The position is not relative to the min/max params.
       * @param name Name of sequence to extract from
       * @param beg beginning index, inclusive
       * @param end ending index, inclusive
       * @return subsequence string
       */
      std::string subseq(const std::string &name,
                         int beg,
                         int end) const {
          int len;
          char *ss = faidx_fetch_seq(_index, name.c_str(), beg, end, &len);
          if (len < 0) {
              throw std::invalid_argument("faidx_fetch_seq error (-2 if c_name not present, -1 general error): " +
                  std::to_string(len));
          }
          std::string ret(ss, len);
          free(ss);
          return ret;
      }

      /**
       * @brief
       * Return sequence length
       * @param name sequence name
       * @return sequence length
       */
      size_t seq_len(std::string name) {
          return faidx_seq_len(_index, name.c_str());
      }

      /**
       * @brief
       * Sequence name given a sequence position.
       * Returns from buffer in WRITE mode.
       * @param i index of sequence
       * @return sequence name
       * @throws std::range_error i is out of sequence index range
       */
      std::string seq_name(const size_t i) const {
          if (i > num_seq()) throw std::range_error("Out of sequence index range.");
          return std::string(faidx_iseq(_index, i));
      }

      /**
       * @brief
       * Get all sequence names in the file or buffer when in write mode.
       * @return vector of sequence names.
       */
      const std::vector<std::string> &sequence_names() const {
          return _seq_names;
      }

      /**
       * @brief
       * Return the full set of sequence names and sequences.
       * If in write mode, return the current buffer.
       * @return vector of std::pair<seq_name, seq>
       * @throws std::invalid_argument No file is loaded
       */
      std::vector<std::pair<std::string, std::string>> sequences() const {
          if (!_index) throw std::invalid_argument("No file loaded.");
          std::vector<std::pair<std::string, std::string>> ret;
          for (size_t i = 0; i < num_seq(); ++i) {
              std::string name = std::string(faidx_iseq(_index, i));
              ret.push_back(std::pair<std::string, std::string>(name, seq(name)));
          }
          return ret;
      }

      /**
       * @return true if FASTA index loaded
       */
      bool good() const {
          return _index != 0;
      }

      /**
       * @brief
       * iterator through FASTA records.
       */
      class iter {
        public:


          /**
           * @brief
           * Create an iterator starting at sequence i.
           * @param in ifasta handle
           * @param i index to start: 0 for begin, seq_len() for end
           */
          iter(ifasta &in, size_t i) : _if(in), _i(i), _end(in.num_seq()) {}

          /**
           * @return true if at same index
           */
          bool operator==(const iter &other) const {
              return _i == other._i;
          }

          /**
           * @return true if not at same index
           */
          bool operator!=(const iter &other) const {
              return !operator==(other);
          }

          /*
           * @return true if index is lower than other
           */
          bool operator<(const iter &other) const {
              return _i < other._i;
          }

          /*
           * @return true if index is greater than other
           */
          bool operator>(const iter &other) const {
              return _i > other._i;
          }

          /*
           * @return true if index is lower than / equal other
           */
          bool operator<=(const iter &other) const {
              return _i <= other._i;
          }

          /*
           * @return true if index is greater than / equal other
           */
          bool operator>=(const iter &other) const {
              return _i >= other._i;
          }

          /**
           * @return iterator to next sequence
           */
          iter &operator++() {
              if (_i < _end) ++_i;
              return *this;
          };

          /**
           * @return pair<seq_name, seq>
           */
          std::pair<std::string, std::string> &operator*() {
              _curr.first = _if._seq_names[_i];
              _curr.second = _if.seq(_curr.first);
              return _curr;
          };

          /**
           * @return *pair<seq_name, seq>
           */
          std::pair<std::string, std::string> *operator->() {
              return &operator*();
          };


        private:
          ifasta &_if;
          size_t _i;
          size_t _end;

          std::pair<std::string, std::string> _curr;
      };

      /**
       * @return iterator to first sequence in FASTA file
       */
      iter begin() {
          return iter(*this, 0);
      }

      /**
       * @return iterator to end of FASTA file
       */
      iter end() {
          return iter(*this, num_seq());
      }

      /**
       * @param seq_name sequence name to start at
       * @return iter to seq_name
       */
      iter begin(std::string seq_name) {
          auto f = std::find(_seq_names.begin(), _seq_names.end(), seq_name);
          if (f == _seq_names.end()) return end();
          return iter(*this, f - _seq_names.begin());
      }

    private:
      friend class iter;
      std::vector<std::string> _seq_names;
      std::string _file_name;
      faidx_t *_index = nullptr;
  };
}

TEST_CASE ("FASTA Reading") {
    using std::endl;
    std::string tmpfa = "tmp_tc.fa";

        SUBCASE("Basic read") {
        {
            std::ofstream fao(tmpfa);
            fao
                << ">x" << endl
                << "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTGGTTCCTGGTGCTATGTGTAACTAGTAATGG" << endl
                << "TAATGGATATGTTGGGCTTTTTTCTTTGATTTATTTGAAGTGACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAA" << endl
                << "TGGAAAGAATACAAGATTTGGAGCCAGACAAATCTGGGTTCAAATCCTCACTTTGCCACATATTAGCCATGTGACTTTGA" << endl
                << "ACAAGTTAGTTAATCTCTCTGAACTTCAGTTTAATTATCTCTAATATGGAGATGATACTACTGACAGCAGAGGTTTGCTG" << endl
                << "TGAAGATTAAATTAGGTGATGCTTGTAAAGCTCAGGGAATAGTGCCTGGCATAGAGGAAAGCCTCTGACAACTGGTAGTT" << endl
                << "ACTGTTATTTACTATGAATCCTCACCTTCCTTGACTTCTTGAAACATTTGGCTATTGACCTCTTTCCTCCTTGAGGCTCT" << endl
                << "TCTGGCTTTTCATTGTCAACACAGTCAACGCTCAATACAAGGGACATTAGGATTGGCAGTAGCTCAGAGATCTCTCTGCT" << endl
                << ">y" << endl
                << "GGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTCAAATCCTGGAGCCAGACAAATCTGGGTTC" << endl;
        }
        Vargas::ifasta fa(tmpfa);

            CHECK(fa.num_seq() == 2);
            REQUIRE(fa.sequence_names().size() == 2);
            CHECK(fa.seq_name(0) == "x");
            CHECK(fa.seq_name(1) == "y");
            CHECK(fa.subseq("x", 0, 3) == "CAAA");
            CHECK(fa.subseq("y", 0, 2) == "GGA");
            CHECK(fa.sequence_names()[0] == "x");
            CHECK(fa.sequence_names()[1] == "y");

    }

        SUBCASE("iterator") {
        {
            std::ofstream o(tmpfa);
            o << ">a\nAAA\nAA\n>b\nCCC\nCC\n>c c\nTTT\nTT\n";
        }

            SUBCASE("Normal iterator") {
            Vargas::ifasta fin(tmpfa);
            auto i = fin.begin();

                CHECK(i->first == "a");
                CHECK(i->second == "AAAAA");
            ++i;

                CHECK(i->first == "b");
                CHECK(i->second == "CCCCC");
            ++i;

                CHECK(i->first == "c");
                CHECK(i->second == "TTTTT");
            ++i;

                CHECK(i == fin.end());
            ++i;
                CHECK(i == fin.end());
        }

            SUBCASE("Resuming iterator") {
            Vargas::ifasta fin(tmpfa);
            {
                auto i = fin.begin("B");
                    CHECK(i == fin.end());
            }

            auto i = fin.begin("b");
                CHECK(i->first == "b");
                CHECK(i->second == "CCCCC");
            ++i;

                CHECK(i->first == "c");
                CHECK(i->second == "TTTTT");
            ++i;

                CHECK(i == fin.end());
            ++i;
                CHECK(i == fin.end());

        }
    }

    remove(tmpfa.c_str());
    remove((tmpfa + ".fai").c_str());
}

TEST_CASE ("FASTA Writing") {
        SUBCASE("open constructor") {
        {
            Vargas::ofasta fa("tmp_tc_wr.fa");
            fa.char_per_line(5);
            fa.write("a", "AAAAA");
            fa.write("b", "TT");
            fa.write("c", "CCCCCCCCCCCC");
        }
        std::ifstream in("tmp_tc_wr.fa");
        std::string line;

        std::getline(in, line);
            CHECK(line == ">a");
        std::getline(in, line);
            CHECK(line == "AAAAA");
        std::getline(in, line);
            CHECK(line == ">b");
        std::getline(in, line);
            CHECK(line == "TT");
        std::getline(in, line);
            CHECK(line == ">c");
        std::getline(in, line);
            CHECK(line == "CCCCC");
        std::getline(in, line);
            CHECK(line == "CCCCC");
        std::getline(in, line);
            CHECK(line == "CC");

        remove("tmp_tc_wr.fa");
    }
}

#endif //VARGAS_FASTA_H
