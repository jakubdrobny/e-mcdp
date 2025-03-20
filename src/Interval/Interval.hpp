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
  Interval(std::string chr_name, long long begin, long long end);

  operator std::string() const;
  bool operator<(const Interval &other) const;
  bool operator==(const Interval &other) const;
};

std::ostream &operator<<(std::ostream &os, const Interval &interval);
std::string interval_vector_to_string(std::vector<Interval> &intervals);

#endif // INTERVAL_H
