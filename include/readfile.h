/**
 * Ravi Gaddipati
 * November 24, 2015
 * rgaddip1@jhu.edu
 *
 * Wrapper for a reads file that loads reads and meta information
 * into a Read structure.
 *
 * readfile.h
 */

#ifndef VARGAS_READFILE_H
#define VARGAS_READFILE_H

#include "readsource.h"

namespace vargas {

class ReadFile: public ReadSource {

 public:
  ReadFile() { }
  ReadFile(std::string file) {
    _read_file.open(file);
    if (!_read_file.good()) throw std::invalid_argument("Invalid reads file \"" + file + "\"");
  }
  ~ReadFile() {
    if (_read_file) _read_file.close();
  }

  /*
   * Updates the stored read.
   * @return True if it was successful.
   */
  bool update_read() override;

  /*
   * Comment lines.
   * @return Line delimited comment lines.
   */
  std::string get_header() const override { return header; }

  /*
   * fast forward reads until this line.
   * @param read Raw sequence to resume from (does not return this read).
   */
  void resume_from(std::string read);


 protected:
  std::string line;
  std::ifstream _read_file;
  std::vector<std::string> _split_meta;

};

}

#endif //VARGAS_READFILE_H
