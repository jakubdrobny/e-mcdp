#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <string>

class Logger {
public:
  enum Level { INFO, DEBUG, ERROR };

  explicit Logger(const std::string &output_path = "");
  ~Logger();

  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;
  
  Logger(Logger&& other) noexcept = default;
  Logger& operator=(Logger&& other) noexcept = default;

  void log(Level level, const std::string &message);
  void info(const std::string &message);
  void debug(const std::string &message);
  void error(const std::string &message);

private:
  bool log_to_file = false;
  std::ofstream file_output;

  std::string get_timestamp();
  std::string level_to_string(Level level);
};

#endif // LOGGER_H
