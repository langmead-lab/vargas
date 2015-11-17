//
// Created by gaddra on 11/16/15.
//


#include "../include/graph.h"
#include "googletest/googletest/include/gtest/gtest.h"

TEST(GraphTest, paramConstructor) {
  struct vmatch::Graph::GraphParams gParam;
  gParam.match = 0;
  gParam.mismatch = 1;
  gParam.gap_open = 2;
  gParam.gap_extension = 3;
  vmatch::Graph g(gParam);
  ASSERT_EQ(0, g.getParamsCopy().match);
  ASSERT_EQ(1, g.getParamsCopy().mismatch);
  ASSERT_EQ(2, g.getParamsCopy().gap_open);
  ASSERT_EQ(3, g.getParamsCopy().gap_extension);
}

TEST(GraphTest, defaultParams) {
  vmatch::Graph g;
  ASSERT_EQ(50000, g.getParamsCopy().maxNodeLen);
  ASSERT_EQ(100, g.getParamsCopy().ingroup);
  ASSERT_EQ("", g.getParamsCopy().region);
  ASSERT_EQ("", g.getParamsCopy().buildfile);
  ASSERT_FALSE(g.getParamsCopy().genComplement);
  ASSERT_FALSE(g.getParamsCopy().maxAF);
  ASSERT_EQ(2, g.getParamsCopy().match);
  ASSERT_EQ(2, g.getParamsCopy().mismatch);
  ASSERT_EQ(3, g.getParamsCopy().gap_open);
  ASSERT_EQ(1, g.getParamsCopy().gap_extension);
}

TEST(GraphTest, setScores) {
  vmatch::Graph g;
  g.setScores(1, 2, 3, 4);
  ASSERT_EQ(1, g.getParamsCopy().match);
  ASSERT_EQ(2, g.getParamsCopy().mismatch);
  ASSERT_EQ(3, g.getParamsCopy().gap_open);
  ASSERT_EQ(4, g.getParamsCopy().gap_extension);
}

TEST(GraphTest, setParams) {
  vmatch::Graph g;
  struct vmatch::Graph::GraphParams gParam = g.getParamsCopy();
  gParam.match = 0;
  gParam.mismatch = 1;
  gParam.gap_open = 2;
  gParam.gap_extension = 3;
  g.setParams(gParam);
  ASSERT_EQ(0, g.getParamsCopy().match);
  ASSERT_EQ(1, g.getParamsCopy().mismatch);
  ASSERT_EQ(2, g.getParamsCopy().gap_open);
  ASSERT_EQ(3, g.getParamsCopy().gap_extension);
}
