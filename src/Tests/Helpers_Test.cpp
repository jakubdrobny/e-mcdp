#include "../Helpers/Helpers.h"
#include <gtest/gtest.h>

TEST(MergeNonDisjointIntervalsTest, EmptyVector) {
  std::vector<Interval> intervals;
  std::vector<Interval> expected;
  EXPECT_EQ(merge_non_disjoint_intervals(intervals), expected);
}

TEST(MergeNonDisjointIntervalsTest, SingleInterval) {
  std::vector<Interval> intervals = {{"chr1", 1, 5}};
  std::vector<Interval> expected = {{"chr1", 1, 5}};
  EXPECT_EQ(merge_non_disjoint_intervals(intervals), expected);
}

TEST(MergeNonDisjointIntervalsTest, MergingOverlappingIntervals) {
  std::vector<Interval> intervals = {{"chr1", 1, 5}, {"chr1", 3, 8}};
  std::vector<Interval> expected = {{"chr1", 1, 8}};
  EXPECT_EQ(merge_non_disjoint_intervals(intervals), expected);
}

TEST(MergeNonDisjointIntervalsTest, NonOverlappingIntervals) {
  std::vector<Interval> intervals = {{"chr1", 1, 5}, {"chr1", 6, 8}};
  std::vector<Interval> expected = {{"chr1", 1, 5}, {"chr1", 6, 8}};
  EXPECT_EQ(merge_non_disjoint_intervals(intervals), expected);
}

TEST(MergeNonDisjointIntervalsTest,
     EqualEndAndBeginShouldMergeBecauseWeNeedGaps) {
  std::vector<Interval> intervals = {{"chr1", 1, 5}, {"chr1", 5, 8}};
  std::vector<Interval> expected = {{"chr1", 1, 8}};
  EXPECT_EQ(merge_non_disjoint_intervals(intervals), expected);
}

TEST(MergeNonDisjointIntervalsTest, MixedChromosomes) {
  std::vector<Interval> intervals = {{"chr1", 1, 5}, {"chr2", 1, 5}};
  std::vector<Interval> expected = {{"chr1", 1, 5}, {"chr2", 1, 5}};
  EXPECT_EQ(merge_non_disjoint_intervals(intervals), expected);
}
