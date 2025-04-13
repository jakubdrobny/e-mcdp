#include "../Args/Args.hpp"
#include "../Interval/Interval.hpp"
#include "../Model/WindowModel.hpp"
#include <gtest/gtest.h>

class WindowModelRunTest : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    ref_intervals = {{"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 5124, 5224},
                     {"chr1", 6333, 6433}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
    query_intervals = {
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
  }

  static std::vector<Interval> ref_intervals;
  static std::vector<Interval> query_intervals;
};

// Initialize static members
std::vector<Interval> WindowModelRunTest::ref_intervals;
std::vector<Interval> WindowModelRunTest::query_intervals;

TEST_F(WindowModelRunTest, NonOverlappingWindows) {
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

TEST_F(WindowModelRunTest, EmptySectionMerge) {
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

TEST_F(WindowModelRunTest, EmptyWindow) {
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

TEST_F(WindowModelRunTest, OverflowingIntervals) {
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

TEST_F(WindowModelRunTest, NestedWindows) {
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

TEST_F(WindowModelRunTest, LongOverflow) {
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

TEST_F(WindowModelRunTest, G24_TEST) {
  std::vector<Interval> ref_ints = {
      {"chr1", 223, 323}, {"chr1", 956, 1056}, {"chr1", 2146, 2246}, {"chr1", 2678, 2778}};
  std::vector<Interval> query_ints = {
      {"chr1", 42, 142},    {"chr1", 297, 397},   {"chr1", 547, 647},   {"chr1", 874, 974},   {"chr1", 1151, 1251},
      {"chr1", 1263, 1363}, {"chr1", 1434, 1534}, {"chr1", 1587, 1687}, {"chr1", 1733, 1833}, {"chr1", 1861, 1961},
      {"chr1", 2079, 2179}, {"chr1", 2325, 2425}, {"chr1", 2486, 2586}, {"chr1", 2678, 2778}, {"chr1", 2899, 2999}};
  std::vector<Interval> windows = {{"chr1", 0, 2000}, {"chr1", 1000, 3000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 3000}};
  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}

TEST_F(WindowModelRunTest, SmallTest) {
  std::vector<Interval> ref_ints = {{"chr1", 956, 1056}};
  std::vector<Interval> query_ints = {
      {"chr1", 42, 142},    {"chr1", 297, 397},   {"chr1", 547, 647},   {"chr1", 874, 974},   {"chr1", 1151, 1251},
      {"chr1", 1263, 1363}, {"chr1", 1434, 1534}, {"chr1", 1587, 1687}, {"chr1", 1733, 1833}, {"chr1", 1861, 1961}};
  std::vector<Interval> windows = {{"chr1", 0, 1000}, {"chr1", 1000, 2000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 2000}};
  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST).run();
  std::vector<WindowResult> resultsTest =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::TEST).run();
  ASSERT_EQ(resultsNaive, resultsFast);
  ASSERT_EQ(resultsNaive, resultsTest);
}

TEST(LargeWindowModelTest, DISABLED_LargeTests) {
  Args args(logger);
  args.ref_intervals_file_path = "data/02-synth-data/g24_8.ref.tsv";
  args.query_intervals_file_path = "data/02-synth-data/g24_8.query.tsv";
  args.chr_size_file_path = "data/02-synth-data/g24_sizes.tsv";
  args.windows_source = "basic";
  args.windows_size = 2000;
  args.windows_step = 1000;

  std::vector<Interval> ref_intervals = load_intervals(args.ref_intervals_file_path);
  std::vector<Interval> query_intervals = load_intervals(args.query_intervals_file_path);
  std::unordered_map<std::string, long long> chr_sizes = load_chr_sizes(args.chr_size_file_path);

  std::unordered_set<std::string> chr_names = load_chr_names_from_chr_sizes(chr_sizes);

  ref_intervals = filter_intervals_by_chr_name(ref_intervals, chr_names);
  query_intervals = filter_intervals_by_chr_name(query_intervals, chr_names);

  ref_intervals = merge_non_disjoint_intervals(ref_intervals);
  query_intervals = merge_non_disjoint_intervals(query_intervals);

  ref_intervals = remove_empty_intervals(ref_intervals);
  query_intervals = remove_empty_intervals(query_intervals);

  for (std::string windows_source : {"basic", "dense"}) {
    args.windows_source = windows_source;
    std::vector<Interval> windows = load_windows(args, chr_sizes);
    windows = filter_intervals_by_chr_name(windows, chr_names);
    windows = remove_empty_intervals(windows);

    std::vector<WindowResult> resultsNaive =
        WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::NAIVE).run();
    std::vector<WindowResult> resultsFast =
        WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::FAST).run();
    std::vector<WindowResult> resultsTest =
        WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::TEST).run();

    for (int i = 0; i < resultsNaive.size(); i++) {
      ASSERT_EQ(resultsNaive[i], resultsFast[i]);
      ASSERT_EQ(resultsFast[i], resultsTest[i]);
    }
  }
}
