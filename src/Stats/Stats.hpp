#ifndef STATS_H
#define STATS_H

#include "../Enums/Enums.hpp"
#include "../Results/WindowResult.hpp"

class Stats {
public:
  Stats(WindowResult result, Significance significance);

  long long get_overlap_count();
  Interval get_window();
  std::vector<long double> get_probs();
  long double get_pvalue();
  long double get_mean();
  long double get_variance();
  long double get_standard_deviation();
  long double get_zscore();
  Significance get_significance();

private:
  long long overlap_count;
  Interval window;
  std::vector<long double> probs;
  long double pvalue, mean, variance, standard_deviation, zscore;
  Significance significance;

  long double calculate_pvalue();
  long double calculate_mean();
  long double calculate_variance();
  long double calculate_standard_deviation();
  long double calculate_zscore();
};

#endif // STATS_H
