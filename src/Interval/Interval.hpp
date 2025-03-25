#ifndef INTERVAL_H
#define INTERVAL_H

#include <ostream>
#include <string>
#include <vector>

// half open intervals [b, e)
class Interval {
public:
  std::string chr_name;
  long long begin, end;

  Interval();
  Interval(const std::string &chr_name, const long long &begin, const long long &end);

  operator std::string() const;
  bool operator<(const Interval &other) const;
  bool operator==(const Interval &other) const;

  long long length();
};

std::ostream &operator<<(std::ostream &os, const Interval &interval);
std::string interval_vector_to_string(std::vector<Interval> &intervals);

#endif // INTERVAL_H
