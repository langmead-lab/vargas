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
    readfile.open(file);
    if (!readfile.good()) throw std::invalid_argument("Invalid read file.");
  }
  ~ReadFile() {
    if (readfile) readfile.close();
  }

  /*
   * Gets the current read.
   * @return Read that has already been generated.
   */
  Read &getRead() override { return read; }

  /*
   * Updates the stored read.
   * @return True if it was successful.
   */
  bool updateRead() override;

  /*
   * Comment lines.
   * @return Line delimited comment lines.
   */
  std::string getHeader() const override { return header; }


 protected:
  std::string line;
  std::ifstream readfile;
  std::vector<std::string> splitMeta;

};

}

#endif //VARGAS_READFILE_H
