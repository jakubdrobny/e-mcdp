#include "Interval.hpp"
#include <string>

Interval::Interval() : chr_name(), begin(), end() {}

Interval::Interval(const std::string &chr_name, const long long &begin, const long long &end)
    : chr_name(chr_name), begin(begin), end(end) {}

Interval::operator std::string() const {
  return chr_name + ": [" + std::to_string(begin) + ", " + std::to_string(end) + ")";
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

long long Interval::length() const { return end - begin; }

std::string Interval::get_chr_name() const { return this->chr_name; }

long long Interval::get_begin() const { return this->begin; }

long long Interval::get_end() const { return this->end; }

std::ostream &operator<<(std::ostream &os, const Interval &interval) {
  os << interval.chr_name + ": [" + std::to_string(interval.begin) + ", " + std::to_string(interval.end) + ")";
  return os;
}

std::string interval_vector_to_string(std::vector<Interval> &intervals) {
  std::string result;
  for (size_t idx = 0; idx < intervals.size(); idx++) {
    if (idx > 0)
      result += ",";
    result += "{" + std::string(intervals[idx]) + "}";
  }
  return result;
}
