#ifndef WINDOWRESULT_H
#define WINDOWRESULT_H

#include "../Interval/Interval.hpp"
#include <array>
#include <vector>

using MultiProbs = std::array<std::array<std::vector<long double>, 2>, 2>;

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

  Interval get_window() const;
  long long get_overlap_count() const;
  std::vector<long double> get_probs() const;
  MultiProbs get_multi_probs() const;

  bool operator==(const WindowResult &other) const;
};

std::ostream &operator<<(std::ostream &os, const WindowResult &result);

#endif // WINDOWRESULT_H
