//
// Created by gaddra on 11/16/15.
//


#include "../include/graph.h"
#include "googletest/googletest/include/gtest/gtest.h"

TEST(graphTest, test1) {
  vmatch::Graph g("data/r5", "data/v5", "out");
  g.exportDOT("out.dot");
}

// Empty class to test protected member
/**
class GraphTestR: public vmatch::Graph, public ::testing::Test { };

//TODO fatal and graphBuild test fixtures are equivalent, how to use the same one for multiple tests?
class fatalGraphBuildTests: public ::testing::TestWithParam<int>, public vmatch::Graph {
 protected:
  int testNo;
  virtual void SetUp() {
    testNo = GetParam();
    ref.open("data/r" + std::to_string(testNo));
    vcf.open("data/v" + std::to_string(testNo));
    cb.open("data/b" + std::to_string(testNo));
    cg.open("data/g" + std::to_string(testNo));
    ASSERT_TRUE(ref.good());
    ASSERT_TRUE(vcf.good());
  }
  virtual void TearDown() {
    ref.close();
    vcf.close();
    cb.close();
    cg.close();
  }
  void build() {
    buildout << buildGraph(ref, vcf).rdbuf();
    correctbuild << cb.rdbuf();

    buildGraph(buildout);
    exportDOT(out);
    correctout << cg.rdbuf();
  }
  std::stringstream buildout, correctbuild;
  std::stringstream out, correctout;
  std::ifstream ref, vcf, cb, cg;
};
class graphBuildTests: public ::testing::TestWithParam<int> {
 protected:
  int testNo;
  virtual void SetUp() {
    testNo = GetParam();
    ref = ("data/r" + std::to_string(testNo));
    vcf = ("data/v" + std::to_string(testNo));
    cb = ("data/b" + std::to_string(testNo));
    cg = ("data/g" + std::to_string(testNo));
  }
  virtual void TearDown() {

  }
  void build() {
    vmatch::Graph g(ref, vcf, "out");
  }
  std::stringstream buildout, correctbuild;
  std::stringstream out, correctout;
  std::string ref, vcf, cb, cg;
};

**/
//TODO exportDOT tests
/**
//********************* Max Node len tests *********************
TEST(GraphTestNodeLen, maxNodelen) {
  struct vmatch::Graph::GraphParams p;
  p.maxNodeLen = 1;
  vmatch::Graph g(p);
  std::ifstream ref, vcf, cb, cg;
  std::stringstream buildout, correctbuild, out, correctout;

  ref.open("data/r8");
  vcf.open("data/v8");
  cb.open("data/b8");
  cg.open("data/g8");
  ASSERT_TRUE(ref.good());
  ASSERT_TRUE(vcf.good());

  buildout << g.buildGraph(ref, vcf).rdbuf();
  correctbuild << cb.rdbuf();
  ASSERT_EQ(correctbuild.str(), buildout.str());

  g.buildGraph(buildout);
  g.exportDOT(out);
  correctout << cg.rdbuf();
  ASSERT_EQ(correctout.str(), out.str());

}

//********************* Graph Output tests *********************

TEST_P(fatalGraphBuildTests, FatalTests) {
  ASSERT_ANY_THROW(build());
}
INSTANTIATE_TEST_CASE_P(fatalGraphBuild, fatalGraphBuildTests, ::testing::Range(1, 2));


TEST_P(graphBuildTests, graphs) {
  build();
  // Check buildfile output
  ASSERT_EQ(correctbuild.str(), buildout.str());
  // Check graph output
  ASSERT_EQ(correctout.str(), out.str());
}
INSTANTIATE_TEST_CASE_P(graphBuild, graphBuildTests, ::testing::Range(2, 7));


//********************* Region tests *********************
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
**/
/********************* Constructor tests *********************/
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

/********************* Params test *********************/
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
  struct vmatch::Graph::GraphParams gParam;
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
