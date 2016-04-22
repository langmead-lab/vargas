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

namespace vargas {

class Xcoder {
 public:
  Xcoder() { simdvbyteinit(); }
  ~Xcoder() { }

  /**
   * data is a list of ints, length is the number of ints, out is the compressed data
   * @param data data to compress
   * @param length length of data
   * @param out pointer to pointer of compressed output
   * @return the length of out in bytes
   *
   */
  size_t compress(const uint32_t *data, size_t length, uint8_t **out);

  /**
   * data is the compressed input, length is the length of data in bytes, and out is the list of ints
   * @param data compressed data
   * @param length length of compressed data
   * @param out pointer to pointer of decompressed data
   * /return the number of ints recovered
   */
  size_t inflate(const uint8_t *data, size_t length, uint32_t **out);

  size_t compress(const std::vector<uint32_t> &vec, uint8_t **out);

  size_t inflate(const uint8_t *data, size_t length, std::vector<uint32_t> &vec);


  /**
   * Encodes the input stream into plaintext characters. Used to convert compressed data to plaintext.
   */
  void encode(std::istream &in, std::ostream &out) {
    E.encode(in, out);
  }

  /**
   * Decodes the plaintext into the ostream.
   */
  void decode(std::istream &in, std::ostream &out) {
    D.decode(in, out);
  }

  std::string encode(const uint8_t *data, size_t len);

  size_t decode(std::string &in, uint8_t **out);

  std::string compressAndEncode(const std::vector<uint32_t> &vec);

  bool testCompression(std::vector<uint32_t> vec);


 protected:
  base64::encoder E;
  base64::decoder D;

};

}

#endif //VARGAS_XCODER_H
