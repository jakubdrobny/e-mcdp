#include "MarkovChain.hpp"
#include "../Logger/Logger.hpp"

#include <iostream>

MarkovChain::MarkovChain(TransitionMatrix T, TransitionMatrix T_MOD) : T(T), T_MOD(T_MOD) {
  this->calculate_stationary_distribution();
}

MarkovChain::MarkovChain(TransitionMatrix T, TransitionMatrix T_MOD, std::array<long double, 2> stationary_distribution)
    : T(T), T_MOD(T_MOD) {}

MarkovChain::MarkovChain(long long chr_len, const std::vector<Interval> &query_intervals) {
  this->calculate_transition_matrices(chr_len, query_intervals);
  this->calculate_stationary_distribution();
}

TransitionMatrix MarkovChain::get_T() const { return this->T; }

TransitionMatrix MarkovChain::get_T_MOD() const { return this->T_MOD; }

StationaryDistribution MarkovChain::get_stationary_distribution() const { return this->stationary_distribution; };

void MarkovChain::print() const {
  std::cout << "T:";
  for (int i : {0, 1})
    std::cout << this->T[i][0] << " " << this->T[i][1] << "\n";
  std::cout << "T_MOD:";
  for (int i : {0, 1})
    std::cout << this->T_MOD[i][0] << " " << this->T_MOD[i][1] << "\n";
  std::cout << "stationary_distribution: " << this->stationary_distribution[0] << " "
            << this->stationary_distribution[1] << "\n";
}

// calculaters T transition matrix
void MarkovChain::calculate_base_transition_matrix(long long chr_size, const std::vector<Interval> &query_intervals) {
  long double L = chr_size;
  long double len_Q = query_intervals.size();
  long double weight_Q = 0;
  for (Interval interval : query_intervals)
    weight_Q += interval.length();

  this->T[0][1] = (len_Q) / (L - weight_Q - 1);
  this->T[0][0] = 1 - this->T[0][1];

  this->T[1][0] = (len_Q) / (weight_Q);
  this->T[1][1] = 1 - this->T[1][0];
}

// calculates T and T_mod (its just T, with zeros in second col)
void MarkovChain::calculate_transition_matrices(long long chr_size, const std::vector<Interval> &query_intervals) {
  if (query_intervals.empty()) {
    logger.error("Query intervals should not be empty.");
    exit(1);
  }

  calculate_base_transition_matrix(chr_size, query_intervals);

  this->T_MOD[0][0] = this->T[0][0];
  this->T_MOD[1][0] = this->T[1][0];
}

void MarkovChain::calculate_stationary_distribution() {
  // b = mat[0][1], a = 1 - b, c = mat[1][0], d = 1 - c
  // ax + cy = x and bx + dy = y and pi_0 + pi_1 = 1 should hold
  // solving for pi_0 and pi_1 we get:
  // we can derive that pi_0 = c/(b+c) and pi_1 = b/(b+c)
  long double b = this->T[0][1], c = this->T[1][0];
  long double denom = b + c;

  // is not irreducible
  if (std::abs(denom) < 1e-9) {
    logger.error("Can't calculate stationary distribution. P[0][1] or P[1][0] are 0.");
    exit(1);
  }

  this->stationary_distribution = {c / denom, b / denom};
}
