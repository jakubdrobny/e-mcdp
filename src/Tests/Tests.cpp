#include "../Helpers/Helpers.hpp"
#include "../Interval/Interval.hpp"
#include "../Model/WindowModel.hpp"
#include <csignal>
#include <gtest/gtest-death-test.h>
#include <gtest/gtest.h>
#include <math.h>

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

TEST(MergeNonDisjointIntervalsTest, EqualEndAndBeginShouldMergeBecauseWeNeedGaps) {
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
  std::array<std::array<long double, 2>, 2> P = {{{{0.7L, 0.3L}}, {{0.2L, 0.8L}}}};
  auto pi = MarkovChain(P, {}).get_stationary_distribution();
  EXPECT_NEAR(pi[0], 0.4L, 1e-10L);
  EXPECT_NEAR(pi[1], 0.6L, 1e-10L);
}

TEST(GetStationaryDistributionTest, AbsorbingState0) {
  std::array<std::array<long double, 2>, 2> P = {{{{1.0L, 0.0L}}, {{0.2L, 0.8L}}}};
  auto pi = MarkovChain(P, {}).get_stationary_distribution();
  EXPECT_NEAR(pi[0], 1.0L, 1e-10L);
  EXPECT_NEAR(pi[1], 0.0L, 1e-10L);
}

TEST(GetStationaryDistributionTest, AbsorbingState1) {
  std::array<std::array<long double, 2>, 2> P = {{{{0.5L, 0.5L}}, {{0.0L, 1.0L}}}};
  auto pi = MarkovChain(P, {}).get_stationary_distribution();
  EXPECT_NEAR(pi[0], 0.0L, 1e-10L);
  EXPECT_NEAR(pi[1], 1.0L, 1e-10L);
}

TEST(GetStationaryDistributionTest, NotIrreducible) {
  std::array<std::array<long double, 2>, 2> P = {{{{1.0L, 0.0L}}, {{0.0L, 1.0L}}}};
  EXPECT_EXIT(MarkovChain(P, {}).get_stationary_distribution(), testing::ExitedWithCode(1), "");
}

TEST(GetStationaryDistributionTest, SumsToOne) {
  std::array<std::array<long double, 2>, 2> P = {{{{0.9L, 0.1L}}, {{0.4L, 0.6L}}}};
  auto pi = MarkovChain(P, {}).get_stationary_distribution();
  long double sum = pi[0] + pi[1];
  EXPECT_NEAR(sum, 1.0L, 1e-10L);
}

TEST(GetWindowsIntervalsTest, NonOverlappingWindows) {
  std::vector<Interval> intervals = {{"", 1, 10}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 15}, {"", 15, 50}};
  auto windows_intervals = WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{{"", 1, 10}}, {{"", 20, 30}}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, OverlappingWindows) {
  std::vector<Interval> intervals = {{"", 10, 12}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 15}, {"", 5, 50}};
  auto windows_intervals = WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{{"", 10, 12}}, {{"", 10, 12}, {"", 20, 30}}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, SlicesIntervalsCorrectly) {
  std::vector<Interval> intervals = {{"", 10, 12}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 11}, {"", 5, 25}};
  auto windows_intervals = WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{{"", 10, 11}}, {{"", 10, 12}, {"", 20, 25}}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, EmptySlicedIntervalsAreExcluded) {
  std::vector<Interval> intervals = {{"", 10, 12}, {"", 20, 30}};
  std::vector<Interval> windows = {{"", 0, 10}, {"", 30, 50}};
  auto windows_intervals = WindowModel::get_windows_intervals(windows, intervals);
  std::vector<std::vector<Interval>> expected = {{}, {}};
  EXPECT_EQ(expected, windows_intervals);
}

TEST(GetWindowsIntervalsTest, FailOnOverlappingIntervals) {
  std::vector<Interval> intervals = {{"", 0, 10}, {"", 5, 15}};
  std::vector<Interval> windows = {{"", 5, 10}, {"", 3, 50}};
  EXPECT_EXIT(WindowModel::get_windows_intervals(windows, intervals), testing::ExitedWithCode(1), "");
}

TEST(GetWindowsIntervalsTest, NoIntervalsInWindow) {
  std::vector<Interval> intervals = {{"", 10, 20}};
  std::vector<Interval> windows = {{"", 20, 30}};
  EXPECT_EQ(WindowModel::get_windows_intervals(windows, intervals), std::vector<std::vector<Interval>>{{}});
}

TEST(GetWindowsIntervalsTest, NestedWindow) {
  std::vector<Interval> intervals = {{"", 7, 12}};
  std::vector<Interval> windows = {{"", 5, 10}, {"", 0, 20}};
  std::vector<std::vector<Interval>> expected = {{{"", 7, 10}}, {{"", 7, 12}}};
  EXPECT_EQ(WindowModel::get_windows_intervals(windows, intervals), expected);
}

TEST(SplitWindowsTest, EmptyInput) {
  std::vector<Interval> windows;
  WindowSectionSplitResult result = split_windows_into_non_overlapping_sections(windows, {}, {});
  EXPECT_TRUE(result.get_sections().empty());
  EXPECT_TRUE(result.get_spans().empty());
}

TEST(SplitWindowsTest, SingleWindow) {
  std::vector<Interval> windows = {{"chr1", 100, 200}};
  WindowSectionSplitResult result = split_windows_into_non_overlapping_sections(windows, {}, {});

  ASSERT_EQ(result.get_sections().size(), 1);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));

  ASSERT_EQ(result.get_spans().size(), 1);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 1));
}

TEST(SplitWindowsTest, NonOverlappingWindows) {
  std::vector<Interval> windows = {{"chr1", 100, 200}, {"chr1", 300, 400}};
  WindowSectionSplitResult result = split_windows_into_non_overlapping_sections(windows, {}, {});

  ASSERT_EQ(result.get_sections().size(), 2);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));
  EXPECT_EQ(result.get_sections()[1], Interval("chr1", 300, 400));

  ASSERT_EQ(result.get_spans().size(), 2);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 1));
  EXPECT_EQ(result.get_spans()[1], Interval("chr1", 1, 2));
}

TEST(SplitWindowsTest, OverlappingWindows) {
  std::vector<Interval> windows = {{"chr1", 100, 300}, {"chr1", 200, 400}};
  WindowSectionSplitResult result = split_windows_into_non_overlapping_sections(windows, {}, {});

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
  WindowSectionSplitResult result = split_windows_into_non_overlapping_sections(windows, {}, {});

  ASSERT_EQ(result.get_sections().size(), 2);
  EXPECT_EQ(result.get_sections()[0], Interval("chr1", 100, 200));
  EXPECT_EQ(result.get_sections()[1], Interval("chr1", 200, 300));

  ASSERT_EQ(result.get_spans().size(), 2);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 1));
  EXPECT_EQ(result.get_spans()[1], Interval("chr1", 1, 2));
}

TEST(SplitWindowsTest, ComplexOverlapping) {
  std::vector<Interval> windows = {{"chr1", 100, 400}, {"chr1", 200, 500}, {"chr1", 300, 600}};
  WindowSectionSplitResult result = split_windows_into_non_overlapping_sections(windows, {}, {});

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

TEST(SplitWindowsTest, IntervalOverflow) {
  std::vector<Interval> windows = {{"chr1", 100, 200}, {"chr1", 200, 300}};
  WindowSectionSplitResult result = split_windows_into_non_overlapping_sections(
      windows, {{"chr1", 100, 250}, {"chr1", 270, 300}}, {{"chr1", 100, 250}, {"chr1", 270, 300}});

  ASSERT_EQ(result.get_sections().size(), 2);
  EXPECT_EQ(result.get_sections()[0], Section("chr1", 100, 200, false, true, false, true));
  EXPECT_EQ(result.get_sections()[1], Section("chr1", 200, 300, true, false, true, false));

  ASSERT_EQ(result.get_spans().size(), 2);
  EXPECT_EQ(result.get_spans()[0], Interval("chr1", 0, 1));
  EXPECT_EQ(result.get_spans()[1], Interval("chr1", 1, 2));
}

std::vector<Interval> ref_intervals = {{"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 5124, 5224},
                                       {"chr1", 6333, 6433}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
std::vector<Interval> query_intervals = {
    {"chr1", 40, 140},    {"chr1", 213, 313},   {"chr1", 505, 605},   {"chr1", 638, 738},   {"chr1", 859, 959},
    {"chr1", 996, 1096},  {"chr1", 1135, 1235}, {"chr1", 1328, 1428}, {"chr1", 1448, 1548}, {"chr1", 1728, 1828},
    {"chr1", 1863, 1963}, {"chr1", 2159, 2259}, {"chr1", 2278, 2378}, {"chr1", 2472, 2572}, {"chr1", 2652, 2752},
    {"chr1", 2783, 2883}, {"chr1", 2935, 3035}, {"chr1", 3102, 3202}, {"chr1", 3207, 3307}, {"chr1", 3311, 3411},
    {"chr1", 3468, 3568}, {"chr1", 3572, 3672}, {"chr1", 3780, 3880}, {"chr1", 3935, 4035}, {"chr1", 4331, 4431},
    {"chr1", 4582, 4682}, {"chr1", 4854, 4954}, {"chr1", 5079, 5179}, {"chr1", 5218, 5318}, {"chr1", 5443, 5543},
    {"chr1", 5650, 5750}, {"chr1", 6050, 6150}, {"chr1", 6198, 6298}, {"chr1", 6590, 6690}, {"chr1", 6739, 6839},
    {"chr1", 6898, 6998}, {"chr1", 7288, 7388}, {"chr1", 7428, 7528}, {"chr1", 7682, 7782}, {"chr1", 8166, 8266},
    {"chr1", 8399, 8499}, {"chr1", 8548, 8648}, {"chr1", 8755, 8855}, {"chr1", 9140, 9240}, {"chr1", 9242, 9342},
    {"chr1", 9388, 9488}, {"chr1", 9647, 9747}, {"chr1", 9749, 9849}, {"chr1", 9914, 10000}};

TEST(WindowModelRunTest, NonOverlappingWindows) {
  std::vector<Interval> windows = {
      {"chr1", 0, 2000}, {"chr1", 2000, 4000}, {"chr1", 4000, 6000}, {"chr1", 6000, 8000}, {"chr1", 8000, 10000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}

TEST(WindowModelRunTest, EmptySectionMerge) {
  std::vector<Interval> windows = {{"chr1", 0, 10000}, {"chr1", 5000, 15000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 15000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}

TEST(WindowModelRunTest, EmptyWindow) {
  std::vector<Interval> windows = {{"chr1", 0, 5000}, {"chr1", 5000, 10000}, {"chr1", 10000, 15000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 15000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}

TEST(WindowModelRunTest, OverflowingIntervals) {
  std::vector<Interval> _intervals = {
      {"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 3000, 7000}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
  std::vector<Interval> windows = {
      {"chr1", 0, 4000}, {"chr1", 2000, 6000}, {"chr1", 4000, 8000}, {"chr1", 6000, 10000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}

TEST(WindowModelRunTest, NestedWindows) {
  std::vector<Interval> _intervals = {
      {"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 3000, 7000}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
  std::vector<Interval> windows = {{"chr1", 2000, 6000}, {"chr1", 2700, 5400}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}

TEST(WindowModelRunTest, LongOverflow) {
  std::vector<Interval> _intervals = {
      {"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 3000, 9000}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
  std::vector<Interval> windows = {{"chr1", 0, 6000}, {"chr1", 2000, 8000}, {"chr1", 4000, 10000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}
