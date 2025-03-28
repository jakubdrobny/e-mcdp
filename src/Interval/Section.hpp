#ifndef SECTION_H
#define SECTION_H

#include "Interval.hpp"
#include <string>

// half open intervals [b, e), but with extra info (whether first interval is from the previous section and wheter last
// interval is going into the next section)
class Section : public Interval {
public:
  bool first_interval_intersected, last_interval_intersected;

  Section(const std::string &chr_name, const long long &begin, const long long &end, bool first_interval_intersected,
          bool last_interval_intersected);
};

#endif // INTERVAL_H
