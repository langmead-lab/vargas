//
// Created by gaddra on 11/17/15.
//

#include "../include/utils.h"
#include "googletest/googletest/include/gtest/gtest.h"

using std::vector;

TEST(splitTest, splitEmpty) {
  vector<std::string> a = split("", ',');
  ASSERT_EQ(0, a.size());
}
TEST(splitTest, splitDelim) {
  vector<std::string> a = split(",", ',');
  ASSERT_EQ(0, a.size());
}
TEST(splitTest, splitSpaceDelim) {
  vector<std::string> a = split(" abcd ", ' ');
  ASSERT_EQ(1, a.size());
  ASSERT_EQ("abcd", a[0]);
}
TEST(splitTest, splitEndDelim) {
  vector<std::string> a = split("abcd,", ',');
  ASSERT_EQ(1, a.size());
  ASSERT_EQ("abcd", a[0]);
}
TEST(splitTest, splitBeginDelim) {
  vector<std::string> a = split(",abcd", ',');
  ASSERT_EQ(1, a.size());
  ASSERT_EQ("abcd", a[0]);
}
TEST(splitTest, splitTwo) {
  vector<std::string> a = split("ab,cd", ',');
  ASSERT_EQ(2, a.size());
  ASSERT_EQ("ab", a[0]);
  ASSERT_EQ("cd", a[1]);
}
TEST(splitTest, splitMultipleDelim) {
  vector<std::string> a = split("ab,,cd", ',');
  ASSERT_EQ(2, a.size());
  ASSERT_EQ("ab", a[0]);
  ASSERT_EQ("cd", a[1]);
}
