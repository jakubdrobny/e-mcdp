#ifndef MODEL_H
#define MODEL_H

#include "../Helpers/Helpers.hpp"
#include "../Interval/Interval.hpp"
#include "../MarkovChain/MarkovChain.hpp"

#include <functional>
#include <vector>

class Model {
public:
  std::vector<Interval> ref_intervals, query_intervals;
  ChrSizesVector chr_sizes;
  std::string method;

  using ProbMethod =
      std::function<std::vector<long double>(std::vector<Interval>, std::vector<Interval>, MarkovChain, long long)>;
  ProbMethod prob_method;

  Model();
  Model(std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals, ChrSizesMap chr_sizes_map,
        std::string method);

  long double eval_pvalue(long long overlap_count);

  static std::vector<long double> eval_probs_single_chr_direct(std::vector<Interval> ref_intervals,
                                                               std::vector<Interval> query_intervals,
                                                               const MarkovChain &markov_chain, long long chr_size);
  static std::vector<long double> eval_probs_single_chr_direct_eigen(std::vector<Interval> ref_intervals,
                                                                     std::vector<Interval> query_intervals,
                                                                     const MarkovChain &markov_chain,
                                                                     long long chr_size);
  static std::vector<std::vector<std::vector<long double>>>
  eval_probs_single_chr_direct_new(const std::vector<Interval> &ref_intervals, long long window_start,
                                   long long window_end, const MarkovChain &markov_chain);

protected:
  static std::vector<Interval> select_intervals_by_chr_name(std::vector<Interval> &intervals, size_t &intervals_idx,
                                                            std::string chr_name);
};

#endif // MODEL_H
