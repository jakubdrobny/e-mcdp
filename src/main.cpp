#include "Args/Args.hpp"
#include "Enums/Enums.hpp"
#include "Helpers/Helpers.hpp"
#include "Logger/Logger.hpp"
#include "Model/Model.hpp"
#include "Model/WindowModel.hpp"
#include "Output/Output.hpp"
#include "Timer/Timer.hpp"

#include <chrono>
#include <gtest/gtest.h>
#include <string>
#include <unordered_map>
#include <vector>

#define RUN_TESTS 1

Logger logger;

int main(int argc, char *argv[]) {
  Timer timer;

  if (RUN_TESTS) {
    ::testing::InitGoogleTest();
    if (RUN_ALL_TESTS()) {
      logger.error("Some tests have failed. Please see the log above.");
      return 1;
    }
    return 0;
  }

  logger = Logger();

  Args args(logger);
  args.parse_args(argc, argv);
  args.debug_args();
  logger.info("Further logs will be in the file specified by the --o flag.");

  if (args.log_file_path != "")
    logger = Logger(args.log_file_path);

  Output output(args.output_file_path);

  logger.info("Loading reference interval set from: " + args.ref_intervals_file_path);
  std::vector<Interval> ref_intervals = load_intervals(args.ref_intervals_file_path);

  logger.info("Loading query interval set from: " + args.query_intervals_file_path);
  std::vector<Interval> query_intervals = load_intervals(args.query_intervals_file_path);

  logger.info("Loading chromosome sizes from: " + args.chr_size_file_path);
  std::unordered_map<std::string, long long> chr_sizes = load_chr_sizes(args.chr_size_file_path);

  size_t raw_ref_count = ref_intervals.size();
  size_t raw_query_count = query_intervals.size();

  std::unordered_set<std::string> chr_names = load_chr_names_from_chr_sizes(chr_sizes);

  ref_intervals = filter_intervals_by_chr_name(ref_intervals, chr_names);
  query_intervals = filter_intervals_by_chr_name(query_intervals, chr_names);

  ref_intervals = merge_non_disjoint_intervals(ref_intervals);
  query_intervals = merge_non_disjoint_intervals(query_intervals);

  ref_intervals = remove_empty_intervals(ref_intervals);
  query_intervals = remove_empty_intervals(query_intervals);

  if (args.statistic == Statistic::BASES) {
    ref_intervals = split_intervals_into_ones(ref_intervals);
    query_intervals = split_intervals_into_ones(query_intervals);
  }

  logger.info("Number of reference intervals: " + std::to_string(ref_intervals.size()) + " (" +
              std::to_string(raw_ref_count) + " before merging)");
  logger.info("Number of query intervals: " + std::to_string(query_intervals.size()) + " (" +
              std::to_string(raw_query_count) + " before merging)");
  logger.info("Number of chromosomes: " + std::to_string(chr_sizes.size()));

  long long overlap_count = count_overlaps(ref_intervals, query_intervals);
  logger.info("Overlap count: " + std::to_string(overlap_count));

  if (!args.windows_source.empty()) {
    // ideme pocitat pre okna
    logger.info("Loading window sizes...");
    std::vector<Interval> windows = load_windows(args, chr_sizes);
    long long raw_window_count = windows.size();
    windows = filter_intervals_by_chr_name(windows, chr_names);
    windows = remove_empty_intervals(windows);

    logger.info("Number of windows: " + std::to_string(windows.size()) + " (" + std::to_string(raw_window_count) +
                " before preprocessing)");

    WindowModel model(windows, ref_intervals, query_intervals, chr_sizes, args.algorithm);
    std::vector<WindowResult> results = model.run();

    output.print("chr_name\tbegin\tend\toverlap_count\tp-value\n");
    for (WindowResult result : results) {
      long double p_value = calculate_joint_pvalue({result.get_probs()}, result.get_overlap_count());
      Interval window = result.get_window();
      output.print(window.chr_name + "\t" + std::to_string(window.begin) + "\t" + std::to_string(window.end) + "\t" +
                   std::to_string(result.get_overlap_count()) + "\t" + std::to_string(p_value) + "\n");
    }

    long double duration = timer.elapsed<std::chrono::milliseconds>();
    logger.debug("Time taken to calculate p-value: " + std::to_string(duration) + " milliseconds\n");
  } else {
    // ideme pocitat pre cely genom spolu
    Model model(ref_intervals, query_intervals, chr_sizes, args.method);

    // measure time from here, previous parts are just tests and i/o
    long double p_value = model.eval_pvalue(overlap_count);
    long double duration = timer.elapsed<std::chrono::milliseconds>();

    logger.info("eval_pvalue: p-value=" + to_string(p_value));

    logger.debug("Time taken to calculate p-value: " + std::to_string(duration) + " milliseconds\n");
  }

  return 0;
}
