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
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlow);
  ASSERT_EQ(resultsNaive, resultsFastBad);
  ASSERT_EQ(resultsNaive, resultsFast);
}

TEST_F(WindowModelRunTest, EmptySectionMerge) {
  std::vector<Interval> windows = {{"chr1", 0, 10000}, {"chr1", 5000, 15000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 15000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlow);
  ASSERT_EQ(resultsNaive, resultsFastBad);
  ASSERT_EQ(resultsNaive, resultsFast);
}

TEST_F(WindowModelRunTest, EmptyWindow) {
  std::vector<Interval> windows = {{"chr1", 0, 5000}, {"chr1", 5000, 10000}, {"chr1", 10000, 15000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 15000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlow);
  ASSERT_EQ(resultsNaive, resultsFastBad);
  ASSERT_EQ(resultsNaive, resultsFast);
}

TEST_F(WindowModelRunTest, OverflowingIntervals) {
  std::vector<Interval> _intervals = {
      {"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 3000, 7000}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
  std::vector<Interval> windows = {
      {"chr1", 0, 4000}, {"chr1", 2000, 6000}, {"chr1", 4000, 8000}, {"chr1", 6000, 10000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlow);
  ASSERT_EQ(resultsNaive, resultsFastBad);
  ASSERT_EQ(resultsNaive, resultsFast);
}

TEST_F(WindowModelRunTest, NestedWindows) {
  std::vector<Interval> _intervals = {
      {"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 3000, 7000}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
  std::vector<Interval> windows = {{"chr1", 2000, 6000}, {"chr1", 2700, 5400}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlow);
  ASSERT_EQ(resultsNaive, resultsFastBad);
  ASSERT_EQ(resultsNaive, resultsFast);
}

TEST_F(WindowModelRunTest, LongOverflow) {
  std::vector<Interval> _intervals = {
      {"chr1", 1909, 2009}, {"chr1", 2694, 2794}, {"chr1", 3000, 9000}, {"chr1", 9299, 9399}, {"chr1", 9629, 9729}};
  std::vector<Interval> windows = {{"chr1", 0, 6000}, {"chr1", 2000, 8000}, {"chr1", 4000, 10000}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};

  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, _intervals, query_intervals, chr_sizes_map, Algorithm::FAST).run();
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlow);
  ASSERT_EQ(resultsNaive, resultsFastBad);
  ASSERT_EQ(resultsNaive, resultsFast);
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
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST).run();
  ASSERT_EQ(resultsNaive, resultsSlowBad);
  ASSERT_EQ(resultsNaive, resultsSlow);
  ASSERT_EQ(resultsNaive, resultsFastBad);
  ASSERT_EQ(resultsNaive, resultsFast);
}

TEST_F(WindowModelRunTest, SmallTest) {
  std::vector<Interval> ref_ints = {
      {"chr1", 4254, 4354}, {"chr1", 4391, 4491}, {"chr1", 4551, 4651}, {"chr1", 4779, 4879}, {"chr1", 5046, 5146},
      {"chr1", 5365, 5465}, {"chr1", 5506, 5606}, {"chr1", 5667, 5767}, {"chr1", 5792, 5892}, {"chr1", 6139, 6239}};
  std::vector<Interval> query_ints = {{"chr1", 4284, 4384}, {"chr1", 4465, 4565}, {"chr1", 4721, 4821},
                                      {"chr1", 4890, 4990}, {"chr1", 5073, 5173}, {"chr1", 5714, 5814},
                                      {"chr1", 6150, 6250}, {"chr1", 6261, 6361}};
  std::vector<Interval> windows = {
      {"chr1", 2200, 4300}, {"chr1", 2250, 4350}, {"chr1", 2300, 4400}, {"chr1", 2350, 4450}, {"chr1", 2400, 4500},
      {"chr1", 2450, 4550}, {"chr1", 2500, 4600}, {"chr1", 2550, 4650}, {"chr1", 2600, 4700}, {"chr1", 2650, 4750},
      {"chr1", 2700, 4800}, {"chr1", 2750, 4850}, {"chr1", 2800, 4900}, {"chr1", 2850, 4950}, {"chr1", 2900, 5000},
      {"chr1", 2950, 5050}, {"chr1", 3000, 5100}, {"chr1", 3050, 5150}, {"chr1", 3100, 5200}, {"chr1", 3150, 5250},
      {"chr1", 3200, 5300}, {"chr1", 3250, 5350}, {"chr1", 3300, 5400}, {"chr1", 3350, 5450}, {"chr1", 3400, 5500},
      {"chr1", 3450, 5550}, {"chr1", 3500, 5600}, {"chr1", 3550, 5650}, {"chr1", 3600, 5700}, {"chr1", 3650, 5750},
      {"chr1", 3700, 5800}, {"chr1", 3750, 5850}, {"chr1", 3800, 5900}, {"chr1", 3850, 5950}, {"chr1", 3900, 6000},
      {"chr1", 3950, 6050}, {"chr1", 4000, 6100}, {"chr1", 4050, 6150}, {"chr1", 4100, 6200}, {"chr1", 4150, 6250},
      {"chr1", 4200, 6300}, {"chr1", 4250, 6350}, {"chr1", 4300, 6400}};
  ChrSizesMap chr_sizes_map = {{"chr1", 10000}};
  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST).run();
  for (size_t i = 0; i < resultsNaive.size(); i++) {
    ASSERT_EQ(resultsNaive[i], resultsSlowBad[i]);
    ASSERT_EQ(resultsNaive[i], resultsSlow[i]);
    ASSERT_EQ(resultsNaive[i], resultsFastBad[i]);
    ASSERT_EQ(resultsNaive[i], resultsFast[i]);
  }
}

TEST_F(WindowModelRunTest, SmallTest2) {
  std::vector<Interval> ref_ints = {{"chr1", 26759, 26859}, {"chr1", 27335, 27435}, {"chr1", 27611, 27711},
                                    {"chr1", 27813, 27913}, {"chr1", 28005, 28105}, {"chr1", 28150, 28250}};
  std::vector<Interval> query_ints = {{"chr1", 26226, 26326}, {"chr1", 26386, 26486}, {"chr1", 26571, 26671},
                                      {"chr1", 26765, 26865}, {"chr1", 26996, 27096}, {"chr1", 27192, 27292},
                                      {"chr1", 27387, 27487}, {"chr1", 27518, 27618}, {"chr1", 27718, 27818},
                                      {"chr1", 27922, 28022}, {"chr1", 28072, 28172}, {"chr1", 28209, 28309}};
  std::vector<Interval> windows = {
      {"chr1", 24050, 26150}, {"chr1", 24100, 26200}, {"chr1", 24150, 26250}, {"chr1", 24200, 26300},
      {"chr1", 24250, 26350}, {"chr1", 24300, 26400}, {"chr1", 24350, 26450}, {"chr1", 24400, 26500},
      {"chr1", 24450, 26550}, {"chr1", 24500, 26600}, {"chr1", 24550, 26650}, {"chr1", 24600, 26700},
      {"chr1", 24650, 26750}, {"chr1", 24700, 26800}, {"chr1", 24750, 26850}, {"chr1", 24800, 26900},
      {"chr1", 24850, 26950}, {"chr1", 24900, 27000}, {"chr1", 24950, 27050}, {"chr1", 25000, 27100},
      {"chr1", 25050, 27150}, {"chr1", 25100, 27200}, {"chr1", 25150, 27250}, {"chr1", 25200, 27300},
      {"chr1", 25250, 27350}, {"chr1", 25300, 27400}, {"chr1", 25350, 27450}, {"chr1", 25400, 27500},
      {"chr1", 25450, 27550}, {"chr1", 25500, 27600}, {"chr1", 25550, 27650}, {"chr1", 25600, 27700},
      {"chr1", 25650, 27750}, {"chr1", 25700, 27800}, {"chr1", 25750, 27850}, {"chr1", 25800, 27900},
      {"chr1", 25850, 27950}, {"chr1", 25900, 28000}, {"chr1", 25950, 28050}, {"chr1", 26000, 28100},
      {"chr1", 26050, 28150}, {"chr1", 26100, 28200}, {"chr1", 26150, 28250}};
  ChrSizesMap chr_sizes_map = {{"chr1", 30000}};
  std::vector<WindowResult> resultsNaive =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::NAIVE).run();
  std::vector<WindowResult> resultsSlowBad =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::SLOW_BAD).run();
  std::vector<WindowResult> resultsSlow =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::SLOW).run();
  std::vector<WindowResult> resultsFastBad =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST_BAD).run();
  std::vector<WindowResult> resultsFast =
      WindowModel(windows, ref_ints, query_ints, chr_sizes_map, Algorithm::FAST).run();
  for (size_t i = 0; i < resultsNaive.size(); i++) {
    ASSERT_EQ(resultsNaive[i], resultsSlowBad[i]);
    ASSERT_EQ(resultsNaive[i], resultsSlow[i]);
    ASSERT_EQ(resultsNaive[i], resultsFastBad[i]);
    ASSERT_EQ(resultsNaive[i], resultsFast[i]);
  }
}

TEST(LargeWindowModelTest, LargeTests) {
  Args args1(logger);
  args1.ref_intervals_file_path = "data/02-synth-data/g24_8.ref.tsv";
  args1.query_intervals_file_path = "data/02-synth-data/g24_8.query.tsv";
  args1.chr_size_file_path = "data/02-synth-data/g24_sizes.tsv";

  Args args2(logger);
  args2.ref_intervals_file_path = "data/04-synthetic-data/ref_1000_1.tsv";
  args2.query_intervals_file_path = "data/04-synthetic-data/query_1000_1.tsv";
  args2.chr_size_file_path = "data/04-synthetic-data/genomeSize.tsv";

  std::vector<std::vector<std::pair<long long, long long>>> window_confs{
      {{2000, 1000}, {10000, 1000}, {200, 100}, {750, 150}, {2100, 140}}, {{1000000, 10000}}};
  std::vector<Args> args_vec = {args1, args2};
  for (size_t args_idx = 0; args_idx < args_vec.size(); args_idx++) {
    Args args = args_vec[args_idx];

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

    for (auto conf : window_confs[args_idx]) {
      args.windows_size = conf.first;
      args.windows_step = conf.second;

      for (std::string windows_source : {"basic", "dense"}) {
        args.windows_source = windows_source;
        std::vector<Interval> windows = load_windows(args, chr_sizes);
        windows = filter_intervals_by_chr_name(windows, chr_names);
        windows = remove_empty_intervals(windows);

        std::vector<WindowResult> resultsNaive =
            WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::NAIVE).run();
        std::vector<WindowResult> resultsSlowBad =
            WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::SLOW_BAD).run();
        std::vector<WindowResult> resultsSlow =
            WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::SLOW).run();
        std::vector<WindowResult> resultsFastBad =
            WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::FAST_BAD).run();
        std::vector<WindowResult> resultsFast =
            WindowModel(windows, ref_intervals, query_intervals, chr_sizes, Algorithm::FAST).run();

        for (size_t i = 0; i < resultsNaive.size(); i++) {
          ASSERT_EQ(resultsNaive[i], resultsSlowBad[i]);
          ASSERT_EQ(resultsNaive[i], resultsSlow[i]);
          ASSERT_EQ(resultsNaive[i], resultsFastBad[i]);
          ASSERT_EQ(resultsNaive[i], resultsFast[i]);
        }
      }
    }
  }
}
