#ifndef MARKOVCHAIN_H
#define MARKOVCHAIN_H

#include "../Interval/Interval.hpp"
#include <array>
#include <vector>

using TransitionMatrix = std::array<std::array<long double, 2>, 2>;
using StationaryDistribution = std::array<long double, 2>;

class MarkovChain {
public:
  MarkovChain();
  MarkovChain(TransitionMatrix T, TransitionMatrix T_MOD, StationaryDistribution stationary_distribution);
  MarkovChain(TransitionMatrix T, TransitionMatrix T_MOD);
  MarkovChain(long long chr_len, const std::vector<Interval> &query_intervals);

  TransitionMatrix get_T() const;
  TransitionMatrix get_T_MOD() const;
  StationaryDistribution get_stationary_distribution() const;

  void print() const;

private:
  TransitionMatrix T{}, T_MOD{};
  StationaryDistribution stationary_distribution{};

  void calculate_base_transition_matrix(long long chr_size, const std::vector<Interval> &query_intervals);
  void calculate_transition_matrices(long long chr_size, const std::vector<Interval> &query_intervals);
  void calculate_transition_matrices();
  void calculate_stationary_distribution();
};

#endif // MARKOVCHAIN_H
