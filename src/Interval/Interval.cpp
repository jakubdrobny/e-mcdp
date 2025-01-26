#include "Interval.h"
#include <string>

Interval::Interval(): chr_name(), begin(), end() {}

Interval::Interval(std::string chr_name, long long begin, long long end): chr_name(chr_name), begin(begin), end(end) {}

Interval::operator std::string() {
  return chr_name + ": [" + std::to_string(begin) + ", " + std::to_string(end) + ")";
}
