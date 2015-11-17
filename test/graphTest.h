//
// Created by gaddra on 11/16/15.
//

#ifndef VMATCH_GRAPHTEST_H
#define VMATCH_GRAPHTEST_H

#include "../include/graph.h"
#include "googletest/googletest/include/gtest/gtest.h"


class GraphTest : public ::testing::Test {
 public:
  virtual void SetUp() {

  }

  vmatch::Graph g;
};

#endif //VMATCH_GRAPHTEST_H
