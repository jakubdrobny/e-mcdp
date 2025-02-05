#ifndef INTERVAL_H
#define INTERVAL_H

#include <ostream>
#include <string>

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
  std::ostream &operator<<(std::ostream &os);
};

#endif // INTERVAL_H
