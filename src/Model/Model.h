#ifndef MODEL_H
#define MODEL_H

#include "../Helpers/Helpers.h"
#include "../Interval/Interval.h"
#include <functional>
#include <vector>

class Model {
public:
  std::vector<Interval> ref_intervals, query_intervals;
  ChrSizesVector chr_sizes;
  std::string method;

  using ProbMethod = std::function<std::vector<long double>(
      std::vector<Interval>, std::vector<Interval>, long long)>;
  ProbMethod prob_method;

  Model();
  Model(std::vector<Interval> ref_intervals,
        std::vector<Interval> query_intervals, ChrSizesMap chr_sizes_map,
        std::string method);

  long double eval_pvalue(long long overlap_count);
  long long eval_sf();

  static std::vector<long double>
  eval_probs_single_chr_direct(std::vector<Interval> ref_intervals,
                               std::vector<Interval> query_intervals,
                               long long chr_size);
  static std::vector<long double>
  eval_probs_single_chr_direct_eigen(std::vector<Interval> ref_intervals,
                                     std::vector<Interval> query_intervals,
                                     long long chr_size);

private:
  static std::vector<Interval>
  select_intervals_by_chr_name(std::vector<Interval> &intervals,
                               size_t &intervals_idx, std::string chr_name);
};

#endif // MODEL_H
