/**
 * Ravi Gaddipati
 * November 27, 2015
 * rgaddip1@jhu.edu
 *
 * vargas::xcoder is used to compress and encode a list of ints,
 * primarily for use with the list of individuals at each variant.
 *
 * xcoder.h
 */


#ifndef VARGAS_XCODER_H
#define VARGAS_XCODER_H

extern "C" {
#include "../VByte/include/varintdecode.h"
#include "../VByte/include/varintencode.h"
}

#include "../libb64/include/encode.h"
#include "../libb64/include/decode.h"
#include <stdexcept>
#include <string.h>
#include <vector>
#include <sstream>

namespace vargas {

class Xcoder {
 public:
  Xcoder() { simdvbyteinit(); }
  ~Xcoder() { }

  /**
   * data is a list of ints, length is the number of ints, out is the compressed data
   * Returns the length of out in bytes
   */
  size_t compress(uint32_t *data, size_t length, uint8_t **out);

  /**
   * data is the compressed input, length is the length of data in bytes, and out is the list of ints
   * Returns the number of ints recovered
   */
  size_t inflate(const uint8_t *data, size_t length, uint32_t **out);

  size_t compress(std::vector<uint32_t> &vec, uint8_t **out);

  size_t inflate(const uint8_t *data, size_t length, std::vector<uint32_t> &vec);

  void encode(std::istream &in, std::ostream &out) {
    E.encode(in, out);
  }
  void decode(std::istream &in, std::ostream &out) {
    D.decode(in, out);
  }

  std::string encode(const uint8_t *data, size_t len);

  size_t decode(std::string &in, uint8_t **out);

  bool testCompression(std::vector<uint32_t> vec);


 protected:
  base64::encoder E;
  base64::decoder D;

};

}

#endif //VARGAS_XCODER_H
