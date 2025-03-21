#include "../Helpers/Helpers.hpp"
#include "../Interval/Interval.hpp"
#include "../Model/WindowModel.hpp"
#include <csignal>
#include <gtest/gtest-death-test.h>
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

TEST(GetStationaryDistributionTest, StandardCase) {
  std::vector<std::vector<long double>> P = {{0.7L, 0.3L}, {0.2L, 0.8L}};
  auto pi = get_stationary_distribution(P);
  EXPECT_NEAR(pi[0], 0.4L, 1e-10L);
  EXPECT_NEAR(pi[1], 0.6L, 1e-10L);
}

TEST(GetStationaryDistributionTest, AbsorbingState0) {
  std::vector<std::vector<long double>> P = {{1.0L, 0.0L}, {0.2L, 0.8L}};
  auto pi = get_stationary_distribution(P);
  EXPECT_NEAR(pi[0], 1.0L, 1e-10L);
  EXPECT_NEAR(pi[1], 0.0L, 1e-10L);
}

TEST(GetStationaryDistributionTest, AbsorbingState1) {
  std::vector<std::vector<long double>> P = {{0.5L, 0.5L}, {0.0L, 1.0L}};
  auto pi = get_stationary_distribution(P);
  EXPECT_NEAR(pi[0], 0.0L, 1e-10L);
  EXPECT_NEAR(pi[1], 1.0L, 1e-10L);
}

TEST(GetStationaryDistributionTest, NotIrreducible) {
  std::vector<std::vector<long double>> P = {{1.0L, 0.0L}, {0.0L, 1.0L}};
  EXPECT_EXIT(get_stationary_distribution(P), testing::ExitedWithCode(1), "");
}

TEST(GetStationaryDistributionTest, InvalidMatrixSize) {
  std::vector<std::vector<long double>> P = {
      {0.6L, 0.3L, 0.1L}, {0.2L, 0.5L, 0.3L}, {0.1L, 0.4L, 0.5L}};
  EXPECT_EXIT(get_stationary_distribution(P), testing::ExitedWithCode(1), "");
}

TEST(GetStationaryDistributionTest, SumsToOne) {
  std::vector<std::vector<long double>> P = {{0.9L, 0.1L}, {0.4L, 0.6L}};
  auto pi = get_stationary_distribution(P);
  long double sum = pi[0] + pi[1];
  EXPECT_NEAR(sum, 1.0L, 1e-10L);
}

TEST(GetWindowsIntervalsTest, NonOverlappingWindows) {
  std::vector<Interval> intervals = {{"", 1, 10}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 15}, {"", 15, 50}};
  auto windows_intervals =
      WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{{"", 1, 10}}, {{"", 20, 30}}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, OverlappingWindows) {
  std::vector<Interval> intervals = {{"", 10, 12}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 15}, {"", 5, 50}};
  auto windows_intervals =
      WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{{"", 10, 12}},
                                                 {{"", 10, 12}, {"", 20, 30}}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, SlicesIntervalsCorrectly) {
  std::vector<Interval> intervals = {{"", 10, 12}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 11}, {"", 5, 25}};
  auto windows_intervals =
      WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{{"", 10, 11}},
                                                 {{"", 10, 12}, {"", 20, 25}}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, EmptySlicedIntervalsAreExcluded) {
  std::vector<Interval> intervals = {{"", 10, 12}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 10}, {"", 30, 50}};
  auto windows_intervals =
      WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{}, {}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, FailOnOverlappingIntervals) {
  std::vector<Interval> intervals = {{"", 0, 10}, {"", 5, 15}};
  std::vector<Interval> windows = {{"", 5, 10}, {"", 3, 50}};
  EXPECT_EXIT(WindowModel::get_windows_intervals(windows, intervals),
              testing::ExitedWithCode(1), "");
}

TEST(GetWindowsIntervalsTest, NoIntervalsInWindow) {
  std::vector<Interval> intervals = {{"", 10, 20}};
  std::vector<Interval> windows = {{"", 20, 30}};
  EXPECT_EQ(WindowModel::get_windows_intervals(windows, intervals),
            std::vector<std::vector<Interval>>{{}});
}

TEST(GetWindowsIntervalsTest, NestedWindow) {
  std::vector<Interval> intervals = {{"", 7, 12}};
  std::vector<Interval> windows = {{"", 5, 10}, {"", 0, 20}};
  std::vector<std::vector<Interval>> expected = {{{"", 7, 10}}, {{"", 7, 12}}};
  EXPECT_EQ(WindowModel::get_windows_intervals(windows, intervals), expected);
}

TEST(SplitWindowsTest, EmptyInput) {
  std::vector<Interval> windows;
  WindowSectionSplitResult result =
      split_windows_into_non_overlapping_sections(windows);
  EXPECT_TRUE(result.get_sections().empty());
  EXPECT_TRUE(result.get_spans().empty());
}

TEST(SplitWindowsTest, SingleWindow) {
  std::vector<Interval> windows = {{"chr1", 100, 200}};
  WindowSectionSplitResult result =
      split_windows_into_non_overlapping_sections(windows);

  ASSERT_EQ(result.get_sections().size(), 1);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));

  ASSERT_EQ(result.get_spans().size(), 1);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 1));
}

TEST(SplitWindowsTest, NonOverlappingWindows) {
  std::vector<Interval> windows = {{"chr1", 100, 200}, {"chr1", 300, 400}};
  WindowSectionSplitResult result =
      split_windows_into_non_overlapping_sections(windows);

  ASSERT_EQ(result.get_sections().size(), 2);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));
  EXPECT_EQ(result.get_sections()[1], Interval("chr1", 300, 400));

  ASSERT_EQ(result.get_spans().size(), 2);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 1));
  EXPECT_EQ(result.get_spans()[1], Interval("chr1", 1, 2));
}

TEST(SplitWindowsTest, OverlappingWindows) {
  std::vector<Interval> windows = {{"chr1", 100, 300}, {"chr1", 200, 400}};
  WindowSectionSplitResult result =
      split_windows_into_non_overlapping_sections(windows);

  ASSERT_EQ(result.get_sections().size(), 3);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));
  EXPECT_EQ(result.get_sections()[1], Interval("chr1", 200, 300));
  EXPECT_EQ(result.get_sections()[2], Interval("chr1", 300, 400));

  ASSERT_EQ(result.get_spans().size(), 2);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 2));
  EXPECT_EQ(result.get_spans()[1], Interval("chr1", 1, 3));
}

TEST(SplitWindowsTest, AdjacentWindows) {
  std::vector<Interval> windows = {{"chr1", 100, 200}, {"chr1", 200, 300}};
  WindowSectionSplitResult result =
      split_windows_into_non_overlapping_sections(windows);

  ASSERT_EQ(result.get_sections().size(), 2);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));
  EXPECT_EQ(result.get_sections()[1], Interval("chr1", 200, 300));

  ASSERT_EQ(result.get_spans().size(), 2);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 1));
  EXPECT_EQ(result.get_spans()[1], Interval("chr1", 1, 2));
}

TEST(SplitWindowsTest, ComplexOverlapping) {
  std::vector<Interval> windows = {
      {"chr1", 100, 400}, {"chr1", 200, 500}, {"chr1", 300, 600}};
  WindowSectionSplitResult result =
      split_windows_into_non_overlapping_sections(windows);

  ASSERT_EQ(result.get_sections().size(), 5);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));
  EXPECT_EQ(result.get_sections()[1], Interval("chr1", 200, 300));
  EXPECT_EQ(result.get_sections()[2], Interval("chr1", 300, 400));
  EXPECT_EQ(result.get_sections()[3], Interval("chr1", 400, 500));
  EXPECT_EQ(result.get_sections()[4], Interval("chr1", 500, 600));

  ASSERT_EQ(result.get_spans().size(), 3);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 3));
  EXPECT_EQ(result.get_spans()[1], Interval("chr1", 1, 4));
  EXPECT_EQ(result.get_spans()[2], Interval("chr1", 2, 5));
}
