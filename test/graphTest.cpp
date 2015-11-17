//
// Created by gaddra on 11/16/15.
//


#include "../include/graph.h"
#include "googletest/googletest/include/gtest/gtest.h"

class GraphTest : public ::testing::Test {
 public:
  virtual void SetUp() {

  }
  int i = 0;
};

TEST_F(GraphTest, aTest) {
ASSERT_EQ(2, i);
}
