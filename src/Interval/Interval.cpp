#include "Interval.h"
#include <ostream>
#include <string>

Interval::Interval() : chr_name(), begin(), end() {}

Interval::Interval(std::string chr_name, long long begin, long long end)
    : chr_name(chr_name), begin(begin), end(end) {}

Interval::operator std::string() const {
  return chr_name + ": [" + std::to_string(begin) + ", " + std::to_string(end) +
         ")";
}

bool Interval::operator<(const Interval &other) const {
  if (chr_name != other.chr_name)
    return chr_name < other.chr_name;
  if (begin != other.begin)
    return begin < other.begin;
  return end < other.end;
}

bool Interval::operator==(const Interval &other) const {
  return chr_name == other.chr_name && begin == other.begin && end == other.end;
}

std::ostream &Interval::operator<<(std::ostream &os) {
  os << chr_name << ": " << begin << "-" << end;
  return os;
}
