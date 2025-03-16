#include "Args.hpp"
#include <vector>

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
    } else if (flag == "--log") {
      if (i + 1 < argc) {
        log_file_path = argv[++i];
        logger.info("Parse --log: " + log_file_path);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--method") {
      if (i + 1 < argc) {
        method = argv[++i];
        logger.info("Parsed --method: " + output_file_path);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--windows.source") {
      if (i + 1 < argc) {
        windows_source = argv[++i];
        logger.info("Parsed --windows.source: " + windows_source);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--windows.path") {
      if (i + 1 < argc) {
        windows_path = argv[++i];
        logger.info("Parsed --windows.path: " + windows_path);
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--windows.size") {
      if (i + 1 < argc) {
        windows_size = std::stoll(argv[++i]);
        logger.info("Parse --windows.size: " + std::to_string(windows_size));
      } else {
        log_failed_to_parse_args(flag);
      }
    } else if (flag == "--windows.step") {
      if (i + 1 < argc) {
        windows_step = std::stoll(argv[++i]);
        logger.info("Parse --windows.step: " + std::to_string(windows_step));
      } else {
        log_failed_to_parse_args(flag);
      }
    } else {
      log_invalid_arg(flag);
    }
  }

  check_required_args();
}

void Args::debug_args() {
  logger.debug("output: " + output_file_path);
  logger.debug("log: " + log_file_path);
  logger.debug("ref_intervals_file_path: " + ref_intervals_file_path);
  logger.debug("query_intervals_file_path: " + query_intervals_file_path);
  logger.debug("chr_size_file_path: " + chr_size_file_path);
  logger.debug("method: " + method);
  logger.debug("windows.source:" + windows_source);
  logger.debug("windows.path: " + windows_path);
  logger.debug("windows.size:" + std::to_string(windows_size));
  logger.debug("windows.step: " + std::to_string(windows_step));
}

void Args::log_failed_to_parse_args(const std::string &flag) {
  logger.error("Argument " + flag + " is missing.");
  exit(1);
}

void Args::log_invalid_arg(const std::string &flag) {
  logger.error("Invalid argument " + flag + ".");
  exit(1);
}

void Args::check_required_args() {
  std::string missing_args;
  if (query_intervals_file_path.empty())
    missing_args += " --q";
  if (ref_intervals_file_path.empty())
    missing_args += " --r";
  if (chr_size_file_path.empty())
    missing_args += " --chs";

  if (!missing_args.empty()) {
    logger.error("Following arguments are missing:" + missing_args + ".");
    exit(1);
  }
}

void Args::check_invalid_args() {
  if (!windows_source.empty() && windows_source != "file" &&
      windows_source != "basic" && windows_source != "dense") {
    logger.error(
        "--windows.source flag can only have values of file, basic or dense.");
    exit(1);
  }

  if (windows_source == "file" && windows_path.empty()) {
    logger.error(
        "--windows.source set to file, but --windows.path was not set.");
    exit(1);
  }

  if (windows_source == "basic" && windows_size <= 0) {
    logger.error("--windows.source set to basic, but --windows.size was not "
                 "set (or was set to <= 0, which is also invalid)");
    exit(1);
  }

  if (windows_source == "dense") {
    std::vector<std::pair<std::string, bool>> conditions = {
        {"--windows.size", (windows_size <= 0)},
        {"--windows.step", (windows_step <= 0)}};

    for (auto entry : conditions) {
      if (entry.second) {
        logger.error("--windows.source set to dense, but " + entry.first +
                     " was not "
                     "set (or was set to <= 0, which is also invalid)");
        exit(1);
      }
    }
  }
}
