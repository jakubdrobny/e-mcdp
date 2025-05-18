#include "Stats.hpp"
#include "../Helpers/Helpers.hpp"

#include <cmath>

Stats::Stats(WindowResult result, Significance significance) {
  this->significance = significance;
  this->window = result.get_window();
  this->probs = result.get_probs();
  this->overlap_count = result.get_overlap_count();
  this->mean = this->calculate_mean();
  // mean and overlap count need to be calculated before p-value so we can choose enrichment or depletion
  this->pvalue = this->calculate_pvalue();
  this->variance = this->calculate_variance();
  this->standard_deviation = this->calculate_standard_deviation();
  this->zscore = this->calculate_zscore();
}

long double Stats::calculate_pvalue() {
  Significance used_significance = this->significance;
  if (used_significance == Significance::COMBINED) {
    if (this->mean < this->overlap_count) {
      used_significance = Significance::ENRICHMENT;
    } else {
      used_significance = Significance::DEPLETION;
    }
  }
  return calculate_joint_pvalue({this->probs}, this->overlap_count, used_significance);
}

long double Stats::calculate_mean() {
  long double mean = 0;
  for (size_t idx = 0; idx < this->probs.size(); idx++) {
    mean += idx * exp(this->probs[idx]);
  }
  return mean;
}

long double Stats::calculate_variance() {
  long double variance = 0;
  for (size_t idx = 0; idx < this->probs.size(); idx++) {
    variance += (idx - this->mean) * (idx - this->mean) * exp(this->probs[idx]);
  }
  return variance;
}

long double Stats::calculate_standard_deviation() { return std::sqrt(this->variance); }

long double Stats::calculate_zscore() { return (this->overlap_count - this->mean) / this->standard_deviation; }

long long Stats::get_overlap_count() { return this->overlap_count; }
Interval Stats::get_window() { return this->window; }
std::vector<long double> Stats::get_probs() { return this->probs; }
long double Stats::get_pvalue() { return this->pvalue; }
long double Stats::get_mean() { return this->mean; }
long double Stats::get_variance() { return this->variance; }
long double Stats::get_standard_deviation() { return this->standard_deviation; }
long double Stats::get_zscore() { return this->zscore; };
Significance Stats::get_significance() { return this->significance; }
