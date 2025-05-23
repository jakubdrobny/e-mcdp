#ifndef ARGS_H
#define ARGS_H

#include "../Enums/Enums.hpp"
#include "../Logger/Logger.hpp"
#include <string>

class Args {
public:
  Args(Logger &logger);

  void parse_args(int argc, char *argv[]);
  void debug_args();

  std::string chr_size_file_path;
  std::string ref_intervals_file_path;
  std::string query_intervals_file_path;
  std::string output_file_path;
  std::string log_file_path;
  Statistic statistic = Statistic::OVERLAPS;
  Algorithm algorithm = Algorithm::NAIVE;
  Significance significance = Significance::ENRICHMENT;
  std::string windows_source;
  std::string windows_path;
  long long windows_size;
  long long windows_step;
  bool run_tests = false;
  bool show_help = false;

private:
  Logger &logger;

  void log_failed_to_parse_args(const std::string &flag);
  void log_invalid_arg(const std::string &flag);
  void check_required_args();
  void check_invalid_args();
};

#endif // ARGS_H
