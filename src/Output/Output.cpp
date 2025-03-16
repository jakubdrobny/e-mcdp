#include "Output.hpp"
#include "../Logger/Logger.hpp"
#include <fstream>
#include <iostream>
#include <string>

Output::Output(const std::string &filename) : filename(filename) {
  if (!filename.empty()) {
    out_to_file = true;
    out = std::ofstream(filename);
    if (!out) {
      logger.error("Failed to open file at the output path.");
      exit(1);
    }
  }
}

Output::~Output() {
  if (out_to_file && out.is_open()) {
    out.close();
  }
}

void Output::print(const std::string &msg) {
  if (out_to_file) {
    out << msg;
  } else {
    std::cout << msg;
  }
}
