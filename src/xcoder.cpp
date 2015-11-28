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


size_t vargas::Xcoder::compress(const uint32_t *data, size_t length, uint8_t **out) {
  if (length <= 0 || *out != NULL) return 0;
  uint8_t *tmp = (uint8_t *) malloc(length * sizeof(uint32_t));
  size_t len = vbyte_encode(data, length, tmp);

  *out = (uint8_t *) malloc(len * sizeof(uint8_t));
  memcpy(*out, tmp, len);

  free(tmp);
  return len;
}


size_t vargas::Xcoder::inflate(const uint8_t *data, size_t length, uint32_t **out) {
  if (length <= 0 || *out != NULL) return 0;
  uint32_t *tmp = (uint32_t *) malloc(length * sizeof(uint32_t));
  size_t len = masked_vbyte_decode_fromcompressedsize(data, tmp, length);

  *out = (uint32_t *) malloc(len * sizeof(uint32_t));
  memcpy(*out, tmp, len * sizeof(uint32_t));

  free(tmp);
  return len;
}


size_t vargas::Xcoder::compress(const std::vector<uint32_t> &vec, uint8_t **out) {
  return compress(vec.data(), vec.size(), out);
}

size_t vargas::Xcoder::inflate(const uint8_t *data, size_t length, std::vector<uint32_t> &vec) {
  uint32_t *out = NULL;
  size_t num = inflate(data, length, &out);

  vec.clear();
  for (int i = 0; i < num; ++i) {
    vec.push_back(out[i]);
  }
  free(out);
  return num;
}

std::string vargas::Xcoder::encode(const uint8_t *data, size_t len) {
  std::stringstream s, o;
  for (int i = 0; i < len; ++i) {
    s.put(data[i]);
  };
  s.seekg(0);
  encode(s, o);
  return o.str();
}

size_t vargas::Xcoder::decode(std::string &in, uint8_t **out) {
  if (out == NULL || *out != NULL) return 0;

  std::stringstream i(in), o;
  decode(i, o);
  uint32_t N = o.str().size();

  *out = (uint8_t *) malloc(N);

  o.seekg(0);
  char c;
  for (int j = 0; j < N; ++j) {
    o.get(c);
    (*out)[j] = uint8_t(c);
  }

  return N;
}


std::string vargas::Xcoder::compressAndEncode(const std::vector<uint32_t> &vec) {
  if (vec.size() == 0) return "";

  uint8_t *compressed = NULL;
  size_t len = compress(vec, &compressed);
  std::string s = encode(compressed, len);

  free(compressed);

  return s;
}



bool vargas::Xcoder::testCompression(std::vector<uint32_t> vec) {
  if (vec.size() == 0) return false;

  std::vector<uint32_t> out;
  uint32_t N = vec.size();

  uint8_t *compressed = NULL;

  size_t cl = compress(vec, &compressed);
  std::cout << std::endl << cl << " bytes to represent ";

  size_t l = inflate(compressed, cl, out);
  std::cout << l << " 32-bit integers." << std::endl;

  std::cout << "Compression ratio: " << float(cl) / (l * sizeof(uint32_t)) << std::endl;

  if (l != N) {
    std::cerr << "Number of ints doesn't match, vec size: " << N << ", Recoevered: " << l << std::endl;
    free(compressed);
    return false;
  }
  for (int j = 0; j < N; ++j) {
    if (out[j] != vec[j]) {
      std::cerr << "Data values don't match." << std::endl;
      free(compressed);
      return false;
    }
  }
  free(compressed);
  return true;
}