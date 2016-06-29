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
#include "utils.h"
#include "doctest.h"
#include "htslib/faidx.h"


namespace Vargas {
/**
 * @brief
* Provides an interface for a FASTA File. An index is built if one does not
* already exist.
*/
  class FASTAFile {

    public:

      /**
       * @enum MODE
       * Set the FASTAFile handle type. The file can be opened in either
       * read or write mode.
       */
      enum class MODE {
          READ, /**< Open FASTAFile in read mode. */
              WRITE /**< Open FASTAFile in write mode. */
      };

      /**
       * @brief
       * Load a given FASTA and its index, create one if none exists.
       * @param file filename
       * @param mode read or write mode, default read
       */
      FASTAFile(std::string file, MODE mode = MODE::READ) : _file_name(file), _mode(mode) { _init(); }

      /**
       * @brief
       * Create a file handle
       * @mode mode read or Write mode, default read.
       */
      FASTAFile(MODE mode = MODE::READ) : _mode(mode) { }

      ~FASTAFile() {
        close();
      }

      /**
       * @brief
       * Close any opened file and flush any outputs.
       */
      void close() {
        if (_index) fai_destroy(_index);
        _index = nullptr;
        _file_name = "";
        if (_mode == MODE::WRITE) flush(_out);
      }

      /**
       * @brief
       * Open a specified FASTA file and make an index.
       * @param file filename
       * @return -1 on index build error, -2 on open error, 0 otherwise
       */
      void open(const std::string &file) {
        _file_name = file;
        _init();
        if (_mode == MODE::WRITE) flush();
      }

      /**
       * @brief
       * Check the number of sequences in the index. In read mode,
       * get the number of buffered sequences.
       * @return number of sequences in the FASTA file
       */
      size_t num_seq() const {
        if (_mode == MODE::WRITE) return _buffer.size();
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
       * @param beg beginning index
       * @param end ending index, inclusive
       * @return subsequence string
       */
      std::string subseq(const std::string &name, int beg, int end) const {
        if (_mode == MODE::WRITE) throw ("subseq() not valid in WRITE mode.");
        int len;
        char *ss = faidx_fetch_seq(_index, name.c_str(), beg, end, &len);
        std::string ret(ss);
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
        if (_mode == MODE::WRITE) throw ("seq_len() not valid in WRITE mode.");
        return faidx_seq_len(_index, name.c_str());
      }

      /**
       * @brief
       * Sequence name given a sequence position.
       * Returns from buffer in WRITE mode.
       * @param i index of sequence
       * @return sequence name
       */
      std::string seq_name(size_t i) const {
        if (_mode == MODE::WRITE) {
          if (i > _buffer.size()) throw std::range_error("Out of buffer index range.");
          return _buffer[i].second;
        }
        if (i > num_seq()) throw std::range_error("Out of sequence index range.");
        return std::string(faidx_iseq(_index, i));
      }

      /**
       * @brief
       * Get all sequence names in the file or buffer when in write mode.
       * @return vector of sequence names.
       */
      std::vector<std::string> sequence_names() const {
        std::vector<std::string> ret;
        if (_mode == MODE::WRITE) {
          for (auto &p : _buffer) {
            ret.push_back(p.first);
          }
          return ret;
        }
        // READ mode
        for (size_t i = 0; i < num_seq(); ++i) {
          ret.push_back(seq_name(i));
        }
        return ret;
      }

      /**
       * @brief
       * Return the full set of sequence names and sequences.
       * If in write mode, return the current buffer.
       * @return vector of std::pair<seq_name, seq>
       */
      std::vector<std::pair<std::string, std::string>> sequences() const {
        if (_mode == MODE::WRITE) return _buffer;
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
      bool good() const { return _index; }

      /**
       * @return FASTA file name
       */
      std::string file() const { return _file_name; }

      /**
       * @brief
       * Add a sequence to the write buffer.
       * @param name sequence name and any meta information.
       * @param sequence sequence string
       */
      void add(const std::string &name, const std::string &sequence) {
        _buffer.push_back(std::pair<std::string, std::string>(name, sequence));
        if (_out.good() && _buffer.size() >= MAX_BUFFER_SIZE) flush(_out);
      }

      /**
       * @brief
       * Flush the buffer to the output file.
       */
      void flush() {
        if (!_out.good()) throw std::invalid_argument("Invalid output file");
        flush(_out);
      }

      /**
       * @brief
       * Max number of characters per FASTA line
       * @param num Max number of chars.
       */
      void char_per_line(size_t num) {
        if (num > 0) _char_per_line = num;
      }


    protected:

      /**
       * @brief
       * Loads a FASTA index. If the index does not exist, one is created.
       * @details
       * With write mode, an index is not created or used. Return true if
       * output file successfully opened.
       * @return -1 on .fai build error, -2 on index load error. 0 on success
       */
      int _init() {
        if (_mode == MODE::WRITE) {
          if (_file_name.length() != 0) _out.open(_file_name);
          if (!_out.good()) throw std::invalid_argument("Error opening output file.");
          return 0;
        }
        // Check if a Fasta index exists. If it doesn't build it.
        if (!file_exists(_file_name + ".fai")) {
          if (fai_build(_file_name.c_str()) != 0) {
            _file_name = "";
            return -1;
          }
        }
        _index = fai_load(_file_name.c_str());
        if (!_index) return -2;
        return 0;
      }

      /**
       * @brief
       * Flush the buffer to the output stream.
       * @param o output stream to flush to
       */
      void flush(std::ostream &o) {
        if (!o.good()) return;
        for (auto &p : _buffer) {
          o << '>' << p.first << '\n';
          size_t pos = 0;
          while (pos < p.second.length()) {
            o << p.second.substr(pos, _char_per_line) << '\n';
            pos += _char_per_line;
          }
        }
        o << std::flush;
        _buffer.clear();
      }

    private:
      std::string _file_name;
      faidx_t *_index = nullptr;
      const MODE _mode;

      std::vector<std::pair<std::string, std::string>> _buffer;
      std::ofstream _out;

      const size_t MAX_BUFFER_SIZE = 1000;
      size_t _char_per_line = 70;
  };
}

TEST_CASE ("FASTA Reading") {
  using std::endl;
  std::string tmpfa = "tmp_tc.fa";
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
    Vargas::FASTAFile fa(tmpfa);

      CHECK(fa.num_seq() == 2);
      REQUIRE(fa.sequence_names().size() == 2);
      CHECK(fa.seq_name(0) == "x");
      CHECK(fa.seq_name(1) == "y");
      CHECK(fa.subseq("x", 0, 3) == "CAAA");
      CHECK(fa.subseq("y", 0, 2) == "GGA");
      CHECK(fa.sequence_names()[0] == "x");
      CHECK(fa.sequence_names()[1] == "y");

  remove(tmpfa.c_str());
  remove((tmpfa + ".fai").c_str());
}

TEST_CASE ("FASTA Writing") {
      SUBCASE("open constructor") {
    {
        Vargas::FASTAFile fa("tmp_tc_wr.fa", Vargas::FASTAFile::MODE::WRITE);
      fa.char_per_line(5);
      fa.add("a", "AAAAA");
      fa.add("b", "TT");
      fa.add("c", "CCCCCCCCCCCC");
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
