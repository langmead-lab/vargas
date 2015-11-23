//
// Created by gaddra on 11/22/15.
//

#include "../include/vcfStream.h"
#include "googletest/googletest/include/gtest/gtest.h"

// Empty class to test protected member
class vcfstreamTest: public vmatch::vcfstream, public ::testing::Test { };

TEST_F(vcfstreamTest, basicVCF) {
  open("data/v3");
  ASSERT_EQ(fields.chrom, 0);
  ASSERT_EQ(fields.pos, 1);
  ASSERT_EQ(fields.id, 2);
  ASSERT_EQ(fields.ref, 3);
  ASSERT_EQ(fields.alt, 4);
  ASSERT_EQ(fields.qual, 5);
  ASSERT_EQ(fields.filter, 6);
  ASSERT_EQ(fields.info, 7);
  ASSERT_EQ(fields.format, 8);
  ASSERT_EQ(fields.indivOffset, 9);
  ASSERT_EQ(fields.numIndivs, 4);
}

TEST(vcfGetTest, getRecordnoInit) {
  vmatch::vcfstream v;
  vmatch::vcfrecord r;
  ASSERT_ANY_THROW(v.getRecord(r));
}

TEST(vcfGetTest, getRecordSingle) {
  vmatch::vcfstream v("data/v3");
  vmatch::vcfrecord r;
  v.getRecord(r);
  ASSERT_EQ(4, r.pos);
  ASSERT_EQ("A", r.ref);
  // ASSERT_EQ("G", r.altIndivs[0]);
  ASSERT_EQ(1, r.altIndivs.size());
}

TEST(vcfGetTest, getRecordNoRec) {
  vmatch::vcfstream v("data/v3");
  vmatch::vcfrecord r;
  v.getRecord(r);
  ASSERT_FALSE(v.getRecord(r));
}

TEST(vcfGetTest, getRecordTwoAlts) {
  vmatch::vcfstream v("data/v4");
  vmatch::vcfrecord r;
  v >> r >> r;

  std::cerr << r;

  ASSERT_EQ(4, r.pos);
  ASSERT_EQ("A", r.ref);
//  ASSERT_EQ("G", r.altIndivs[0]);
//  ASSERT_EQ("C", r.altIndivs[1]);
  ASSERT_EQ(2, r.altIndivs.size());

}
