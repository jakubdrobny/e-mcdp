#include "Args/Args.h"
#include "Helpers/Helpers.h"
#include "Logger/Logger.h"

#include <gtest/gtest.h>
#include <string>
#include <unordered_map>
#include <vector>

Logger logger;

int main(int argc, char *argv[]) {
  logger = Logger();

  ::testing::InitGoogleTest();
  if (RUN_ALL_TESTS()) {
    logger.error("Some tests have failed. Please see the log above.");
    return 1;
  }

  Args args(logger);
  args.parse_args(argc, argv);
  args.debug_args();
  logger.info("Further logs will be in the file specified by the --o flag.");

  if (args.output_file_path != "")
    logger = Logger(args.output_file_path);

  logger.info("Loading reference interval set from: " +
              args.ref_intervals_file_path);
  std::vector<Interval> ref_intervals =
      load_intervals(args.ref_intervals_file_path);

  logger.info("Loading query interval set from: " +
              args.query_intervals_file_path);
  std::vector<Interval> query_intervals =
      load_intervals(args.query_intervals_file_path);

  logger.info("Loading chromosome sizes from: " + args.chr_size_file_path);
  std::unordered_map<std::string, long long> chr_sizes =
      load_chr_sizes(args.chr_size_file_path);

  size_t raw_ref_count = ref_intervals.size();
  size_t raw_query_count = query_intervals.size();

  std::unordered_set<std::string> chr_names =
      load_chr_names_from_chr_sizes(chr_sizes);

  ref_intervals = filter_intervals_by_chr_name(ref_intervals, chr_names);
  query_intervals = filter_intervals_by_chr_name(query_intervals, chr_names);

  ref_intervals = merge_non_disjoint_intervals(ref_intervals);
  query_intervals = merge_non_disjoint_intervals(query_intervals);

  ref_intervals = remove_empty_intervals(ref_intervals);
  query_intervals = remove_empty_intervals(query_intervals);

  logger.info(
      "Number of reference intervals: " + std::to_string(ref_intervals.size()) +
      " (" + std::to_string(raw_ref_count) + " before merging)");
  logger.info(
      "Number of query intervals: " + std::to_string(query_intervals.size()) +
      " (" + std::to_string(raw_query_count) + " before merging)");
  logger.info("Number of chromosomes: " + std::to_string(chr_sizes.size()));

  long long overlap_count = count_overlaps(ref_intervals, query_intervals);
  logger.info("Overlap count: " + std::to_string(overlap_count));

  return 0;
}
