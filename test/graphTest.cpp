//
// Created by gaddra on 11/16/15.
//


#include "../include/graph.h"
#include "googletest/googletest/include/gtest/gtest.h"
#include "equivilence.h"

// Empty class to test protected member
class GraphTestR: public vargas::Graph, public ::testing::Test { };

// Test fixture to test outputs to given correct output files.
class GraphBuildTests: public ::testing::TestWithParam<int>, public vargas::Graph {
 protected:
  int testNo;
  virtual void SetUp() {
    testNo = GetParam();
    ref.open("data/r" + std::to_string(testNo));
    vcf.open("data/v" + std::to_string(testNo));
    cb.open("data/b" + std::to_string(testNo));
    cg.open("data/g" + std::to_string(testNo));
    ASSERT_TRUE(ref.good());
    ASSERT_TRUE(cb.good());
    ASSERT_TRUE(cg.good());
  }
  virtual void TearDown() {
    ref.close();
    cb.close();
    cg.close();
  }
  void buildAndTest() {
    exportBuildfile(&ref, vcf, buildout);

    buildout.clear();
    buildout.seekg(0, buildout.beg);

    buildGraph(buildout);
    exportDOT(graphout);

    buildout.clear();
    graphout.clear();
    buildout.seekg(0, buildout.beg);
    graphout.seekg(0, graphout.beg);

    ASSERT_TRUE(fileEquiv(cb, buildout, true));
    ASSERT_TRUE(fileEquiv(cg, graphout, false));
  }
  std::stringstream buildout;
  std::stringstream graphout;
  vargas::vcfstream vcf;
  std::ifstream ref, cb, cg;
};

TEST_P(GraphBuildTests, sampleGraph) {
  buildAndTest();
}
INSTANTIATE_TEST_CASE_P(graphBuild, GraphBuildTests, ::testing::Range(2, 8));


//TODO exportDOT tests

/********************* Region tests *********************/
TEST_F(GraphTestR, regionParseEmpty) {
  uint32_t min, max;

  parseRegion("", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(UINT32_MAX, max);
}
TEST_F(GraphTestR, regionParseInvalid) {
  uint32_t min, max;

  parseRegion("ds", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(UINT32_MAX, max);
}
TEST_F(GraphTestR, regionParseNoNum) {
  uint32_t min, max;

  parseRegion(":", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(UINT32_MAX, max);
}
TEST_F(GraphTestR, regionParseNoMin) {
  uint32_t min, max;

  parseRegion(":100", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(UINT32_MAX, max);
}
TEST_F(GraphTestR, regionParseNoMax) {
  uint32_t min, max;

  parseRegion("100:", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(UINT32_MAX, max);
}
TEST_F(GraphTestR, regionParseNegMin) {
  uint32_t min, max;

  parseRegion("-100:", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(UINT32_MAX, max);
}
TEST_F(GraphTestR, regionParseNegMax) {
  uint32_t min, max;

  parseRegion(":-100", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(UINT32_MAX, max);
}
TEST_F(GraphTestR, regionParseCorrect) {
  uint32_t min, max;

  parseRegion("0:5000", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(5000, max);
}
TEST_F(GraphTestR, regionParseDoubleSpace) {
  uint32_t min, max;

  parseRegion("0 : 5000", &min, &max);
  ASSERT_EQ(0, min);
  ASSERT_EQ(5000, max);
}
TEST_F(GraphTestR, regionParseSingleSpace) {
  uint32_t min, max;
  parseRegion("25: 500", &min, &max);
  ASSERT_EQ(25, min);
  ASSERT_EQ(500, max);
}
TEST_F(GraphTestR, regionParseSpaceBegin) {
  uint32_t min, max;
  parseRegion(" 25:500", &min, &max);
  ASSERT_EQ(25, min);
  ASSERT_EQ(500, max);
}

/********************* Constructor tests *********************/
TEST(GraphTest, paramConstructor) {
  struct vargas::Graph::GraphParams gParam;
  gParam.match = 0;
  gParam.mismatch = 1;
  gParam.gap_open = 2;
  gParam.gap_extension = 3;
  vargas::Graph g(gParam);
  ASSERT_EQ(0, g.getParamsCopy().match);
  ASSERT_EQ(1, g.getParamsCopy().mismatch);
  ASSERT_EQ(2, g.getParamsCopy().gap_open);
  ASSERT_EQ(3, g.getParamsCopy().gap_extension);
}

/********************* Params test *********************/
TEST(GraphTest, defaultParams) {
  vargas::Graph g;
  ASSERT_EQ(50000, g.getParamsCopy().maxNodeLen);
  ASSERT_EQ(100, g.getParamsCopy().ingroup);
  ASSERT_EQ("", g.getParamsCopy().region);
  ASSERT_EQ("", g.getParamsCopy().complementSource);
  ASSERT_FALSE(g.getParamsCopy().genComplement);
  ASSERT_FALSE(g.getParamsCopy().maxAF);
  ASSERT_EQ(2, g.getParamsCopy().match);
  ASSERT_EQ(2, g.getParamsCopy().mismatch);
  ASSERT_EQ(3, g.getParamsCopy().gap_open);
  ASSERT_EQ(1, g.getParamsCopy().gap_extension);
}

TEST(GraphTest, setScores) {
  vargas::Graph g;
  g.setScores(1, 2, 3, 4);
  ASSERT_EQ(1, g.getParamsCopy().match);
  ASSERT_EQ(2, g.getParamsCopy().mismatch);
  ASSERT_EQ(3, g.getParamsCopy().gap_open);
  ASSERT_EQ(4, g.getParamsCopy().gap_extension);
}

TEST(GraphTest, setParams) {
  vargas::Graph g;
  struct vargas::Graph::GraphParams gParam;
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
