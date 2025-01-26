#ifndef INTERVAL_H 
#define INTERVAL_H

#include <string>

// half open intervals [b, e)
class Interval {
public:
  std::string chr_name;
  long long begin, end;
  
  Interval();
  Interval(std::string chr_name, long long begin, long long end);

  operator std::string();
};

#endif // INTERVAL_H
