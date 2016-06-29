/**
 * @author Ravi Gaddipati
 * @date November 24, 2015
 * rgaddip1@jhu.edu
 *
 * @brief
 * Wrapper for a reads file that loads reads and meta information
 * into a Read structure.
 *
 * @file
 */

#ifndef VARGAS_READFILE_H
#define VARGAS_READFILE_H

#include "readsource.h"

namespace Vargas {

  /**
   * @brief
   * Interface to a FASTA read file.
   */
class ReadFile: public ReadSource {

 public:
    /**
     * @brief
     * If no file is specified, use stdin.
     */
  ReadFile() { }

    /**
     * @param file filename of reads file
     */
  ReadFile(std::string file) {
    _read_file.open(file);
    if (!_read_file.good()) throw std::invalid_argument("Invalid reads file \"" + file + "\"");
  }

  ~ReadFile() {
    if (_read_file.is_open()) _read_file.close();
  }

    /*
     * @brief
     * Updates the stored read.
     * @return True if it was successful. False if end of file.
     */
  bool update_read() override;

    /*
     * @brief
     * Get comment lines.
     * @return Line delimited comment lines.
     */
  std::string get_header() const override { return header; }

    /*
     * @brief
     * fast forward reads until this line.
     * @param read Raw sequence to resume from (does not return this read).
     */
  void resume_from(std::string read);


 protected:
    std::string line;
    /**< Current line loaded from file */
  std::ifstream _read_file;
  std::vector<std::string> _split_meta;

};

}

#endif //VARGAS_READFILE_H
