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

Logger logger;

int main(int argc, char *argv[]) {
  Timer timer;
  logger = Logger();

  Args args(logger);
  args.parse_args(argc, argv);

  if (args.show_help) {
    logger.info("The program provides a set of flags to operate it:");
    logger.info("--r <path-to-your-ref-intervals-file>\t\t\t- REQUIRED, tells the program where to find the file with "
                "reference annotation intervals");
    logger.info("--q <path-to-your-query-intervals-file>\t\t- REQUIRED, tells the program where to find the file with "
                "query annotation intervals");
    logger.info(
        "--chs <path-to-your-chromosome-sizes-file>\t\t- REQUIRED, tells the program where to find the file with "
        "chromosome sizes");
    logger.info("--log <name-of-log-file>\t\t\t\t- logs will be written into this file, which will be located in the "
                "`data/logs/` directory");
    logger.info(
        "--o <name-of-results-file>\t\t\t\t- results of the program will be written into this file, which will be "
        "located in the `data/output/` directory");
    logger.info("--windows.source <file|basic|dense>\t\t\t- specifies the windows source of windows set");
    logger.info("--windows.path <path-to-your-windows-file>\t\t- required with the `--windows.source file` flag, tells "
                "the program the location of the window set file");
    logger.info(
        "--windows.size <windows-size>\t\t\t\t- required with the `--windows.source <basic|dense>` flags, tells the "
        "program the size of windows to generate");
    logger.info(
        "--windows.step <windows-step>\t\t\t\t- required with the `--windows.source dense` flag, tells the program "
        "the shift when generating overlapping set of windows");
    logger.info("--algorithm <naive|slow_bad|slow|fast_bad|fast>\t- defaults to naive, is used to choose algorithm "
                "when evaluating windows");
    logger.info("--test\t\t\t\t\t\t- if this flag is specified, all other flags (except `--help`) are ignored and all "
                "the tests "
                "in the `src/Tests` are ran and then the program quits");
    logger.info(
        "--help\t\t\t\t\t\t- if this flag is specified, all other flags are ignored and a help text will be shown");
    return 0;
  }

  args.debug_args();

  if (args.run_tests) {
    ::testing::InitGoogleTest();
    if (RUN_ALL_TESTS()) {
      logger.error("Some tests have failed. Please see the log above.");
      return 1;
    }
    return 0;
  }

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
    Model model(ref_intervals, query_intervals, chr_sizes);

    // measure time from here, previous parts are just tests and i/o
    long double p_value = model.eval_pvalue(overlap_count);
    long double duration = timer.elapsed<std::chrono::milliseconds>();

    logger.info("eval_pvalue: p-value=" + to_string(p_value));

    logger.debug("Time taken to calculate p-value: " + std::to_string(duration) + " milliseconds\n");
  }

  return 0;
}
