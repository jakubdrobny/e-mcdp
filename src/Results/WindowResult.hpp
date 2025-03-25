#ifndef WINDOWRESULT_H
#define WINDOWRESULT_H

#include "../Interval/Interval.hpp"
#include <vector>

using MultiProbs = std::vector<std::vector<std::vector<long double>>>;

class WindowResult {
private:
  Interval window;
  long long overlap_count;
  std::vector<long double> probs;
  MultiProbs multi_probs;

public:
  WindowResult();
  WindowResult(Interval window, long long overlap_count, std::vector<long double> probs);
  WindowResult(Interval window, long long overlap_count, MultiProbs multi_probs);

  Interval get_window();
  long long get_overlap_count();
  std::vector<long double> get_probs();
  MultiProbs get_multi_probs();
};

#endif // WINDOWRESULT_H
