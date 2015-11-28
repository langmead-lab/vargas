//
// Created by gaddra on 11/27/15.
//

#include "googletest/googletest/include/gtest/gtest.h"
#include "../include/xcoder.h"

TEST(xcoder, basic) {
  srand(time(NULL));

  std::vector<uint32_t> dat, dat2;
  uint32_t N = 1000;

  for (int i = 0; i < N; ++i) {
    dat.push_back(rand() & 100000);
  }

  vargas::Xcoder x;

  std::string ec = x.compressAndEncode(dat);
  std::cout << std::endl << ec << std::endl;

  uint8_t *out = NULL;
  size_t decodeLen = x.decode(ec, &out);
  size_t lenInflated = x.inflate(out, decodeLen, dat2);

  ASSERT_EQ(N, lenInflated);
  for (int j = 0; j < N; ++j) {
    ASSERT_EQ(dat[j], dat2[j]);
  }

}

