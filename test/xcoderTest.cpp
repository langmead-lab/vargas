//
// Created by gaddra on 11/27/15.
//

#include "googletest/googletest/include/gtest/gtest.h"
#include "../include/xcoder.h"

TEST(xcoder, basic) {

  uint32_t *data, *out = NULL;
  data = new uint32_t[10];

  for (unsigned int i = 0; i < 10; ++i) {
    data[i] = 100;
    std::cout << data[i] << ",";
  }
  std::cout << std::endl;

  uint8_t *compressed = NULL;

  vargas::Xcoder x;

  size_t cl = x.compress(data, 10, &compressed);
  std::cout << cl << std::endl;

  size_t l = x.inflate(compressed, cl, &out);
  std::cout << l << std::endl;

  for (int j = 0; j < 10; ++j) {
    std::cout << uint32_t(out[j]) << ",";
  }

}

