#include "Args.h"

Args::Args(Logger &logger) : logger(logger) {}

void Args::parse_args(int argc, char *argv[]) {
  for (int i = 1; i < argc; i++) {
    std::string flag = argv[i];
    if (flag == "--chs") {
      if (i + 1 < argc) {
        chr_size_file_path = argv[++i];
        logger.info("Parsed --chs: " + chr_size_file_path);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--q") {
      if (i + 1 < argc) {
        query_intervals_file_path = argv[++i];
        logger.info("Parsed --q: " + query_intervals_file_path);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--r") {
      if (i + 1 < argc) {
        ref_intervals_file_path = argv[++i];
        logger.info("Parsed --r: " + ref_intervals_file_path);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--o") {
      if (i + 1 < argc) {
        output_file_path = argv[++i];
        logger.info("Parsed --o: " + output_file_path);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else {
      log_invalid_arg(flag);
    }
  }

  check_missing_args();
}

void Args::debug_args() {
  logger.debug("output: " + output_file_path);
  logger.debug("ref_intervals_file_path: " + ref_intervals_file_path);
  logger.debug("query_intervals_file_path: " + query_intervals_file_path);
  logger.debug("chr_size_file_path: " + chr_size_file_path);
}

void Args::log_failed_to_parse_args(const std::string &flag) {
  logger.error("Argument " + flag + " is missing. Exiting.");
  exit(1);
}

void Args::log_invalid_arg(const std::string &flag) {
  logger.error("Invalid argument " + flag + ". Exiting.");
  exit(1);
}

void Args::check_missing_args() {
  std::string missing_args;
  if (output_file_path.empty())
    missing_args += " --o";
  if (query_intervals_file_path.empty())
    missing_args += " --q";
  if (ref_intervals_file_path.empty())
    missing_args += " --r";
  if (chr_size_file_path.empty())
    missing_args += " --chs";

  if (!missing_args.empty()) {
    logger.error("Following arguments are missing:" + missing_args + ". Exiting.");
    exit(1);
  }
}
