//
// Created by gaddra on 11/27/15.
//

#include "googletest/googletest/include/gtest/gtest.h"
#include "../include/xcoder.h"

TEST(xcoder, basic) {
  srand(time(NULL));

  std::vector<uint32_t> dat, dat2;
  uint32_t N = 500;

  for (int i = 0; i < N; ++i) {
    dat.push_back(rand() & 10000);
  }

  vargas::Xcoder x;

  uint8_t *compressed = NULL, *out = NULL;
  size_t compressedLen = x.compress(dat, &compressed);
  std::string ec = x.encode(compressed, compressedLen);

  size_t decodeLen = x.decode(ec, &out);
  size_t lenInflated = x.inflate(out, decodeLen, dat2);

  std::cout << "Compression ratio: " << float(compressedLen) / (N * sizeof(uint32_t)) << std::endl;

  ASSERT_EQ(N, lenInflated);
  for (int j = 0; j < N; ++j) {
    ASSERT_EQ(dat[j], dat2[j]);
  }

}

