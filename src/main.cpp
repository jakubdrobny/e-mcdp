#include "Args/Args.h"
#include "Logger/Logger.h"
#include "Helpers/Helpers.h"

#include <string>
#include <unordered_map>
#include <vector>

Logger logger;

int main(int argc, char *argv[]) {
  logger = Logger();
  
  Args args(logger);
  args.parse_args(argc, argv);
  args.debug_args();
  logger.info("Further logs will be in the file specified by the --o flag.");

  if (args.output_file_path != "")
    logger = Logger(args.output_file_path);

  logger.info("Loading reference interval set from: " + args.ref_intervals_file_path);
  std::vector<Interval> ref_intervals = load_intervals(args.ref_intervals_file_path);
 
  logger.info("Loading query interval set from: " + args.query_intervals_file_path);
  std::vector<Interval> query_intevals = load_intervals(args.query_intervals_file_path);

  logger.info("Loading chromosome sizes from: " + args.chr_size_file_path);
  std::unordered_map<std::string, long long> chr_sizes = load_chr_sizes(args.chr_size_file_path);

  return 0;
}
