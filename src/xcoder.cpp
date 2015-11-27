/**
 * Ravi Gaddipati
 * November 27, 2015
 * rgaddip1@jhu.edu
 *
 * vargas::xcoder is used to compress and encode a list of ints,
 * primarily for use with the list of individuals at each variant.
 *
 * xcoder.cpp
 */

#include "../include/xcoder.h"

size_t vargas::Xcoder::compress(uint32_t *data, size_t length, uint8_t **out) {
  uint8_t *tmp = new uint8_t[length];
  size_t len = vbyte_encode(data, length, tmp);
  *out = new uint8_t[len];
  memcpy(*out, tmp, len);
  delete[] tmp;
  return len;
}
size_t vargas::Xcoder::inflate(const uint8_t *data, size_t length, uint32_t **out) {
  uint32_t *tmp = new uint32_t[length];
  size_t len = masked_vbyte_decode(data, tmp, length);
  *out = new uint32_t[len];
  memcpy(*out, tmp, len);
  delete[] tmp;
  return len;
}