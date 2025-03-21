#ifndef WINDOWRESULT_H
#define WINDOWRESULT_H

#include "../Interval/Interval.hpp"
#include <vector>

class WindowResult {
private:
  Interval window;
  long long overlap_count;
  std::vector<long double> probs;

public:
  WindowResult();
  WindowResult(Interval window, long long overlap_count,
               std::vector<long double> probs);

  Interval get_window();
  long long get_overlap_count();
  std::vector<long double> get_probs();
};

#endif // WINDOWRESULT_H
