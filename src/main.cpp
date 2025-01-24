#include <iostream>
#include <string>

class Logger {
  enum Level { INFO, DEBUG };

  // leave empty for stdout
  Logger(const std::string &output_path = "") {}
};

class Args {
public:
  std::string chr_size_file_path;        // --chs flag
  std::string ref_intervals_file_path;   // --r flag
  std::string query_intervals_file_path; // --q flag
  std::string output_file_path;          // -o flag

  void parse_args(int argc, char *argv[]) {
    for (int i = 1; i < argc; i++) {
      std::string flag = argv[i];
      if (flag == "--chs") {
        if (i + 1 < argc) {
          chr_size_file_path = argv[i + 1];
        } else {
          log_failed_to_parse_args(flag);
        }
      } else if (flag == "--q") {
        if (i + 1 < argc) {
          query_intervals_file_path = argv[i + 1];
        } else {
          log_failed_to_parse_args(flag);
        }
      } else if (flag == "--r") {
        if (i + 1 < argc) {
          ref_intervals_file_path = argv[i + 1];
        } else {
          log_failed_to_parse_args(flag);
        }
      } else if (flag == "--o") {
        if (i + 1 < argc) {
          output_file_path = argv[i + 1];
        } else {
          log_failed_to_parse_args(flag);
        }
      } else {
        log_invalid_arg(flag);
      }
    }

    check_missing_args();
  }

private:
  void log_failed_to_parse_args(std::string flag) {
    std::cout << "Argument " << flag << " is missing. Exiting\n";
    exit(1);
  }

  void log_invalid_arg(std::string flag) {
    std::cout << "Invalid argument " << flag << ". Exiting\n.";
    exit(1);
  }

  void check_missing_args() {
    std::cout << "Following arguments are missing:";
    if (output_file_path == "")
      std::cout << " --o";
    if (query_intervals_file_path == "")
      std::cout << " --q";
    if (ref_intervals_file_path == "")
      std::cout << " --r";
    if (chr_size_file_path == "")
      std::cout << " --chs";
    std::cout << ".\nExiting.\n";
    exit(1);
  }
};

int main(int argc, char *argv[]) {
  Args args;
  args.parse_args(argc, argv);
  std::cout << "output: " << args.output_file_path << "\n";
  std::cout << "ref_intervals_file_path: " << args.ref_intervals_file_path
            << "\n";
  std::cout << "query_intervals_file_path: " << args.query_intervals_file_path
            << "\n";
  std::cout << "chr_size_file_path: " << args.chr_size_file_path << "\n";

  return 0;
}
