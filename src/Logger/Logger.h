#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <string>

class Logger {
public:
  enum Level { INFO, DEBUG, ERROR };

  explicit Logger(const std::string &output_path = "");
  ~Logger();

  void log(Level level, const std::string &message);

private:
  bool log_to_file = false;
  std::ofstream file_output;

  std::string get_timestamp();
  std::string level_to_string(Level level);
};

#endif // LOGGER_H
