#ifndef SECTION_H
#define SECTION_H

#include "Interval.hpp"
#include <string>

// half open intervals [b, e), but with extra info (whether first interval is from the previous section and wheter last
// interval is going into the next section)
class Section : Interval {
public:
  std::string chr_name;
  long long begin, end;
  bool first_interval_intersected, last_interval_intersected;

  Section(const std::string &chr_name, const long long &begin, const long long &end, bool first_interval_intersected,
          bool last_interval_intersected);

  long long length();
};

#endif // INTERVAL_H
