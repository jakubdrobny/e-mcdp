#include "Logger.h"
#include <ctime>
#include <iostream>
#include <sstream>

Logger::Logger(const std::string &output_path) {
  if (!output_path.empty()) {
    file_output.open(output_path, std::ios::app);
    if (!file_output.is_open()) {
      throw std::runtime_error("Failed to open log file: " + output_path);
    }
    log_to_file = true;
  }
}

Logger::~Logger() {
  if (file_output.is_open()) {
    file_output.close();
  }
}

void Logger::log(Level level, const std::string &message) {
  std::string log_entry =
      "[" + get_timestamp() + "] [" + level_to_string(level) + "] " + message;

  if (log_to_file) {
    file_output << log_entry << std::endl;
  } else {
    std::cout << log_entry << std::endl;
  }
}

std::string Logger::get_timestamp() {
  std::time_t now = std::time(nullptr);
  std::tm *now_tm = std::localtime(&now);

  std::ostringstream oss;
  oss << (now_tm->tm_year + 1900) << "-" << (now_tm->tm_mon + 1) << "-"
      << now_tm->tm_mday << " " << now_tm->tm_hour << ":" << now_tm->tm_min
      << ":" << now_tm->tm_sec;
  return oss.str();
}

std::string Logger::level_to_string(Level level) {
  switch (level) {
  case INFO:
    return "INFO";
  case DEBUG:
    return "DEBUG";
  case ERROR:
    return "ERROR";
  default:
    return "UNKNOWN";
  }
}
