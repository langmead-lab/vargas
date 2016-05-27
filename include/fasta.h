/**
 * @author Ravi Gaddipati (rgaddip1@jhu.edu)
 * @date May 26, 2016
 *
 * Provides a C++ wrapper for htslib. FASTAFile provides an interface
 * to a FASTA formatted file. An index is created for the opened
 * file.
 *
 * @file fasta.h
 */

#ifndef VARGAS_FASTA_H
#define VARGAS_FASTA_H

#include <string>
#include <fstream>
#include "../doctest/doctest/doctest.h"
#include <../htslib/htslib/faidx.h>


namespace vargas {
/**
* Provides an interface for a FASTA File. An index is built if one does not
* already exist.
*/
class FASTAFile {

 public:
  /**
   * Load a given FASTA and index, create if none exists.
   * @param file filename
   */
  FASTAFile(std::string file) : _file_name(file) { _init(); }

  FASTAFile() { }

  ~FASTAFile() {
    close();
  }

  void close() {
    if (_index) fai_destroy(_index);
    _file_name = "";
  }
  /**
   * Open a specified FASTA file and make an index.
   * @param file filename
   * @return -1 on index build error, -2 on open error, 0 otherwise
   */
  int open(std::string file) {
    _file_name = file;
    return _init();
  }

  /**
   * Check the number of sequences in the index.
   * @return number of sequences in the FASTA file
   */
  int num_seq() const { return faidx_nseq(_index); }

  /**
   * Return a subsequence of the FASTA file, absolute 0 based indexing.
   * The position is not relative to the min/max params.
   * @param beg beginning index
   * @param end ending index, inclusive
   * @return subsequence string
   */
  std::string subseq(std::string chr, int beg, int end) const {
    int len;
    char *ss = faidx_fetch_seq(_index, chr.c_str(), beg, end, &len);
    std::string ret(ss);
    free(ss);
    return ret;
  }

  /**
   * @param sequence name
   * @return sequence length
   */
  int seq_len(std::string seq) { return faidx_seq_len(_index, seq.c_str()); }

  /**
   * Sequence name given a sequence ID.
   * @param ID of sequence
   * @return sequence name
   */
  std::string seq_name(int i) const {
    if (i > num_seq()) return "";
    return std::string(faidx_iseq(_index, i));
  }

  /**
   * Get all sequences in the File
   * @return vector of sequence names.
   */
  std::vector<std::string> sequences() const {
    std::vector<std::string> ret;
    for (int i = 0; i < num_seq(); ++i) {
      ret.push_back(seq_name(i));
    }
    return ret;
  }

  bool good() const { return _index != nullptr; }

  std::string file() const { return _file_name; }


 protected:

  /**
   * Loads a FASTA index. If the index does not exist, one is created.
   */
  int _init() {
    // Check if a Fasta index exists. If it doesn't build it.
    bool exists;
    {
      std::ifstream fidx(_file_name + ".fai");
      exists = fidx.good();
    }

    if (!exists) {
      if (fai_build(_file_name.c_str()) < 0) {
        _file_name = "";
        return -1;
      }
    }

    _index = fai_load(_file_name.c_str());
    if (!_index) return -2;
    return 0;
  }

 private:
  std::string _file_name;
  faidx_t *_index = nullptr;
};
TEST_CASE ("FASTA Handler") {
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
  FASTAFile fa(tmpfa);

      CHECK(fa.num_seq() == 2);
      REQUIRE(fa.sequences().size() == 2);
      CHECK(fa.seq_name(0) == "x");
      CHECK(fa.seq_name(1) == "y");
      CHECK(fa.subseq("x", 0, 3) == "CAAA");
      CHECK(fa.subseq("y", 0, 2) == "GGA");
      CHECK(fa.sequences()[0] == "x");
      CHECK(fa.sequences()[1] == "y");

  remove(tmpfa.c_str());
  remove((tmpfa + ".fai").c_str());
}
}

#endif //VARGAS_FASTA_H
