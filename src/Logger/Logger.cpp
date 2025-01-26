#include "Logger.h"
#include <ctime>
#include <iostream>
#include <sstream>
#include <format>

Logger::Logger(const std::string &output_path) {
  if (!output_path.empty()) {
    file_output.open(output_path, std::ios::trunc);
    if (!file_output.is_open()) {
      std::cerr << "Failed to open log file: " + output_path << ". Exiting.\n";
      exit(1);
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
      "[" + get_timestamp() + "] - [" + level_to_string(level) + "]\t- " + message;

  if (log_to_file) {
    file_output << log_entry << std::endl;
  } else {
    std::cout << log_entry << std::endl;
  }
}

void Logger::info(const std::string &message) {
  log(Logger::INFO, message);
}

void Logger::debug(const std::string &message) {
  log(Logger::DEBUG, message);
}

void Logger::error(const std::string &message) {
  log(Logger::ERROR, message);
}

std::string Logger::get_timestamp() {
  std::time_t now = std::time(nullptr);
  std::tm *now_tm = std::localtime(&now);

  std::ostringstream oss;
  oss << (now_tm->tm_year + 1900) << "-" <<
    std::format("{:02}", now_tm->tm_mon + 1) << "-" <<
    std::format("{:02}", now_tm->tm_mday) << " " <<
    std::format("{:02}", now_tm->tm_hour) << ":" <<
    std::format("{:02}", now_tm->tm_min) << ":" <<
    std::format("{:02}", now_tm->tm_sec);
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
