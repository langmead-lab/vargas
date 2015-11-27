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

namespace vargas {

class Xcoder {
 public:
  Xcoder() { simdvbyteinit(); }
  ~Xcoder() { }

  size_t compress(uint32_t *data, size_t length, uint8_t **out);
  size_t inflate(const uint8_t *data, size_t length, uint32_t **out);

  void encode(std::istream &in, std::ostream &out) {
    E.encode(in, out);
  }
  void decode(std::istream &in, std::ostream &out) {
    D.decode(in, out);
  }


 protected:
  base64::encoder E;
  base64::decoder D;

};

}

#endif //VARGAS_XCODER_H
