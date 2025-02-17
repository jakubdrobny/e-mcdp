#include "Model.h"
#include "../Helpers/Helpers.h"
#include "../Logger/Logger.h"
#include <limits>

Model::Model() : ref_intervals(), query_intervals(), chr_sizes() {}

Model::Model(std::vector<Interval> ref_intervals,
             std::vector<Interval> query_intervals, ChrSizesMap chr_sizes_map,
             std::string method)
    : ref_intervals(ref_intervals), query_intervals(query_intervals),
      method(method) {
  sort(ref_intervals.begin(), ref_intervals.end());
  sort(query_intervals.begin(), query_intervals.end());
  chr_sizes = chr_sizes_map_to_array(chr_sizes_map);
  sort(chr_sizes.begin(), chr_sizes.end());

  if (method == "direct") {
    prob_method = Model::eval_probs_single_chr_direct;
  } else if (method == "eigen") {
    prob_method = Model::eval_probs_single_chr_direct_eigen;
  } else {
    logger.error("Invalid method. Given: " + method +
                 ". Expected: [direct, eigen].");
    exit(1);
  }
}

std::vector<Interval>
Model::select_intervals_by_chr_name(std::vector<Interval> &intervals,
                                    size_t &intervals_idx,
                                    std::string chr_name) {
  std::vector<Interval> chr_intervals;
  for (; intervals_idx < intervals.size() &&
         intervals[intervals_idx].chr_name <= chr_name;
       intervals_idx++)
    if (intervals[intervals_idx].chr_name == chr_name)
      chr_intervals.push_back(intervals[intervals_idx]);
  return chr_intervals;
}

long double Model::eval_pvalue(long long overlap_count) {
  std::vector<std::vector<long double>> probs_by_chromosome;

  size_t ref_idx = 0, query_idx = 0;
  for (size_t chr_sizes_idx = 0; chr_sizes_idx < chr_sizes.size();
       chr_sizes_idx++) {
    std::string chr_name = chr_sizes[chr_sizes_idx].first;
    long long chr_size = chr_sizes[chr_sizes_idx].second;

    std::vector<Interval> chr_ref_intervals =
                              Model::select_intervals_by_chr_name(
                                  ref_intervals, ref_idx, chr_name),
                          chr_query_intervals =
                              Model::select_intervals_by_chr_name(
                                  query_intervals, query_idx, chr_name);
    std::vector<long double> probs =
        prob_method(chr_ref_intervals, chr_query_intervals, chr_size);
    probs_by_chromosome.push_back(probs);
  }

  long double joint_pvalue =
      calculate_joint_pvalue(probs_by_chromosome, overlap_count);
  return joint_pvalue;
}

std::vector<long double>
Model::eval_probs_single_chr_direct(std::vector<Interval> ref_intervals,
                                    std::vector<Interval> query_intervals,
                                    long long chr_size) {
  if (ref_intervals.empty() || query_intervals.empty())
    return {0.};

  auto transition_matrics = get_transition_matrices(chr_size, query_intervals);
  nc::NdArray<long double> T = transition_matrics.first,
                           D = transition_matrics.second;

  int m = ref_intervals.size();
  if (ref_intervals[0].begin == 0) {
    logger.warn("First reference interval starts with zero, changing to one!");
    ref_intervals[0].begin = 1;
    if (ref_intervals[0].end - ref_intervals[0].begin == 0) {
      logger.warn("First reference interval has length 0, removing it!");
      ref_intervals.erase(ref_intervals.begin());
    }
  }

  std::vector<Interval> ref_intervals_augmented;
  ref_intervals_augmented.push_back(Interval(
      ref_intervals[0].chr_name, std::numeric_limits<long long>::min(), 0));
  extend(ref_intervals_augmented, ref_intervals);
  ref_intervals_augmented.push_back(
      Interval(ref_intervals[0].chr_name, chr_size,
               std::numeric_limits<long long>::max()));

  nc::NdArray<long double> prev_line(m + 1, 2);
  prev_line(0, 0) = 1;
  nc::NdArray<long double> last_col(m + 1, 2);

  // calculate zero-th row in separate way
  for (int j = 1; j <= m; j++) {
    long long gap =
        ref_intervals_augmented[j].begin - ref_intervals_augmented[j - 1].end;
    if (j == 1)
      gap--;
    if (gap < 0) {
      logger.error("Gap should be non-negative.");
      exit(1);
    }

    long long len =
        ref_intervals_augmented[j].end - ref_intervals_augmented[j].begin;
    if (len < 0) {
      logger.error("Interval length should be non-negative.");
      exit(1);
    }

    nc::NdArray<long double> row_vector = prev_line.row(j);
    row_vector = row_vector.reshape(1, row_vector.size());
    nc::NdArray<long double> result =
        matrix_multiply(row_vector, binary_exponentiation(T, gap));
    for (nc::uint32 k = 0; k < prev_line.numCols(); ++k)
      prev_line(j, k) = result(0, k);
  }

  std::vector<long double> probs;
  return probs;
}

std::vector<long double>
Model::eval_probs_single_chr_direct_eigen(std::vector<Interval> ref_intervals,
                                          std::vector<Interval> query_intervals,
                                          long long chr_size) {
  return {0.};
}
