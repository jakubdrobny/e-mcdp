#include "Model.hpp"
#include "../Logger/Logger.hpp"
#include "../MarkovChain/MarkovChain.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <omp.h>

Model::Model() : ref_intervals(), query_intervals(), chr_sizes() {}

Model::Model(std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals, ChrSizesMap chr_sizes_map,
             std::string method)
    : ref_intervals(ref_intervals), query_intervals(query_intervals), method(method) {
  sort(ref_intervals.begin(), ref_intervals.end());
  sort(query_intervals.begin(), query_intervals.end());
  chr_sizes = chr_sizes_map_to_array(chr_sizes_map);
  sort(chr_sizes.begin(), chr_sizes.end());

  if (method == "direct") {
    prob_method = Model::eval_probs_single_chr_direct;
  } else if (method == "eigen") {
    prob_method = Model::eval_probs_single_chr_direct_eigen;
  } else {
    logger.error("Invalid method. Given: " + method + ". Expected: [direct, eigen].");
    exit(1);
  }
}

std::vector<Interval> Model::select_intervals_by_chr_name(std::vector<Interval> &intervals, size_t &intervals_idx,
                                                          std::string chr_name) {
  std::vector<Interval> chr_intervals;
  for (; intervals_idx < intervals.size() && intervals[intervals_idx].chr_name <= chr_name; intervals_idx++)
    if (intervals[intervals_idx].chr_name == chr_name)
      chr_intervals.push_back(intervals[intervals_idx]);
  return chr_intervals;
}

long double Model::eval_pvalue(long long overlap_count) {
  std::vector<std::vector<long double>> probs_by_chr(chr_sizes.size());
  std::vector<std::vector<Interval>> ref_intervals_by_chr(chr_sizes.size()), query_intervals_by_chr(chr_sizes.size());

  size_t ref_idx = 0, query_idx = 0;
  for (size_t chr_sizes_idx = 0; chr_sizes_idx < chr_sizes.size(); chr_sizes_idx++) {
    std::string chr_name = chr_sizes[chr_sizes_idx].first;

    ref_intervals_by_chr[chr_sizes_idx] = Model::select_intervals_by_chr_name(ref_intervals, ref_idx, chr_name),
    query_intervals_by_chr[chr_sizes_idx] = Model::select_intervals_by_chr_name(query_intervals, query_idx, chr_name);
  }

// sometimes turned off for debugging
#pragma omp parallel for
  for (size_t chr_sizes_idx = 0; chr_sizes_idx < chr_sizes.size(); chr_sizes_idx++) {
    long long chr_size = chr_sizes[chr_sizes_idx].second;
    MarkovChain markov_chain(chr_size, query_intervals_by_chr[chr_sizes_idx]);
    std::vector<long double> probs =
        prob_method(ref_intervals_by_chr[chr_sizes_idx], query_intervals_by_chr[chr_sizes_idx], markov_chain, chr_size);
    probs_by_chr[chr_sizes_idx] = probs;
  }

  long double joint_pvalue = calculate_joint_pvalue(probs_by_chr, overlap_count);
  return joint_pvalue;
}

std::vector<long double> Model::eval_probs_single_chr_direct(std::vector<Interval> ref_intervals,
                                                             std::vector<Interval> query_intervals,
                                                             const MarkovChain &markov_chain, long long chr_size) {
  if (ref_intervals.empty() || query_intervals.empty())
    return {0.};

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
  ref_intervals_augmented.push_back(Interval(ref_intervals[0].chr_name, std::numeric_limits<long long>::min(), 0));
  extend(ref_intervals_augmented, ref_intervals);

  std::vector<std::array<long double, 2>> prev_line(m + 1, std::array<long double, 2>()),
      last_col(m + 1, std::array<long double, 2>());
  prev_line[0][0] = 1;
  StationaryDistribution stationary_distribution = markov_chain.get_stationary_distribution();
  // prev_line[0][0] = stationary_distribution[0];
  // prev_line[0][1] = stationary_distribution[1];

  // calculate zero-th row in separate way
  for (int j = 1; j <= m; j++) {
    long long gap = ref_intervals_augmented[j].begin - ref_intervals_augmented[j - 1].end;
    if (j == 1)
      gap--;
    if (gap < 0) {
      logger.error("Gap should be non-negative.");
      exit(1);
    }

    long long len = ref_intervals_augmented[j].end - ref_intervals_augmented[j].begin;
    if (len < 0) {
      logger.error("Interval length should be non-negative.");
      exit(1);
    }

    std::array<std::array<long double, 2>, 2> result =
        matrix_multiply({{{{prev_line[j - 1][0], prev_line[j - 1][1]}}}},
                        matrix_multiply(binary_exponentiation(markov_chain.get_T(), gap),
                                        binary_exponentiation(markov_chain.get_T_MOD(), len)));
    prev_line[j] = result[0];
  }

  last_col[0] = prev_line[prev_line.size() - 1];

  std::vector<std::array<long double, 2>> next_line(m + 1, std::array<long double, 2>());
  for (int k = 1; k <= m; k++) {
    next_line[k - 1] = {0, 0};

    if (k % 10 == 0) {
      logger.debug("Processing " + std::to_string(k) + "-th line out of " + std::to_string(m) + " rows of DP table...");
    }

    for (int j = k; j < m + 1; j++) {
      long long gap = ref_intervals_augmented[j].begin - ref_intervals_augmented[j - 1].end;
      if (j == 1)
        gap--;
      if (gap < 0) {
        logger.error("Gap should be non-negative.");
        exit(1);
      }

      long long len = ref_intervals_augmented[j].end - ref_intervals_augmented[j].begin;
      if (len < 0) {
        logger.error("Interval length should be non-negative.");
        exit(1);
      }

      // dont_hit = P[j-1, k] * T^gap * D^len
      std::array<std::array<long double, 2>, 2> dont_hit =
          matrix_multiply({{{{next_line[j - 1][0], next_line[j - 1][1]}}}},
                          matrix_multiply(binary_exponentiation(markov_chain.get_T(), gap),
                                          binary_exponentiation(markov_chain.get_T_MOD(), len)));
      // hit = P[j-1, k-1] * T^gap * (T^len - D^l)
      std::array<std::array<long double, 2>, 2> hit =
          matrix_multiply({{{{prev_line[j - 1][0], prev_line[j - 1][1]}}}},
                          matrix_multiply(binary_exponentiation(markov_chain.get_T(), gap),
                                          subtract_matrices(binary_exponentiation(markov_chain.get_T(), len),
                                                            binary_exponentiation(markov_chain.get_T_MOD(), len))));

      // P[j,k] = dont_hit + hit
      next_line[j] = add_matrices(dont_hit, hit)[0];
    }

    last_col[k][0] = next_line[next_line.size() - 1][0];
    last_col[k][1] = next_line[next_line.size() - 1][1];

    prev_line = next_line;
  }

  std::vector<long double> probs(m + 1);
  for (int k = 0; k <= m; k++) {
    long long trailing_gap = chr_size - ref_intervals_augmented[m].end;
    // std::cout << "lastcol: " << last_col[k][0] << ", " << last_col[k][1]
    //<< "\n";
    last_col[k] = matrix_multiply({{last_col[k]}}, binary_exponentiation(markov_chain.get_T(), trailing_gap))[0];
    std::cout << "lastcol: " << last_col[k][0] << ", " << last_col[k][1] << "\n";
    probs[k] = log(last_col[k][0] + last_col[k][1]);
  }

  return probs;
}

template <class T> using SectionProbs = std::vector<std::vector<T>>;

std::vector<std::vector<std::vector<long double>>>
Model::eval_probs_single_chr_direct_new(std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals,
                                        long long chr_size, long long window_start, long long window_end,
                                        const MarkovChain &mc) {
  if (ref_intervals.empty() || query_intervals.empty())
    return {{{0.}, {0.}}, {{0.}, {0.}}};
  std::cout << "chr_size+query_intervals: " << chr_size;
  for (auto x : query_intervals)
    std::cout << " " << x;
  std::cout << "\n";
  int m = ref_intervals.size();
  std::vector<Interval> ref_intervals_augmented;
  ref_intervals_augmented.push_back(
      Interval(ref_intervals[0].chr_name, std::numeric_limits<long long>::min(), window_start));
  extend(ref_intervals_augmented, ref_intervals);
  ref_intervals_augmented.push_back(
      Interval(ref_intervals[0].chr_name, window_end + 1, std::numeric_limits<long long>::max()));

  SectionProbs<std::vector<long double>> probs(2, std::vector<std::vector<long double>>(2));
  for (int start_state : {0, 1}) {
    std::vector<std::vector<long double>> prev_line(m + 1, std::vector<long double>(2)),
        last_col(m + 1, std::vector<long double>(2));
    prev_line[0][start_state] = 1;

    // calculate zero-th row in separate way
    for (int j = 1; j <= m; j++) {
      long long gap = ref_intervals_augmented[j].begin - ref_intervals_augmented[j - 1].end;
      if (gap < 0) {
        logger.error("Gap should be non-negative.");
        exit(1);
      }

      long long len = ref_intervals_augmented[j].end - ref_intervals_augmented[j].begin;
      if (len < 0) {
        logger.error("Interval length should be non-negative.");
        exit(1);
      }

      std::vector<std::vector<long double>> result =
          matrix_multiply(vector_to_2d_matrix(prev_line[j - 1]),
                          matrix_multiply(binary_exponentiation(T, gap), binary_exponentiation(D, len)));
      prev_line[j] = result[0];
    }

    last_col[0] = prev_line[prev_line.size() - 1];

    std::vector<std::vector<long double>> next_line(m + 1, std::vector<long double>(2));
    for (int k = 1; k <= m; k++) {
      next_line[k - 1] = {0, 0};

      if (k % 10 == 0) {
        logger.debug("Processing " + std::to_string(k) + "-th line out of " + std::to_string(m) +
                     " rows of DP table...");
      }

      for (int j = k; j <= m; j++) {
        long long gap = ref_intervals_augmented[j].begin - ref_intervals_augmented[j - 1].end;
        if (gap < 0) {
          logger.error("Gap should be non-negative.");
          exit(1);
        }

        long long len = ref_intervals_augmented[j].end - ref_intervals_augmented[j].begin;
        if (len < 0) {
          logger.error("Interval length should be non-negative.");
          exit(1);
        }

        // dont_hit = P[j-1, k] * T^gap * D^len
        std::vector<std::vector<long double>> dont_hit =
            matrix_multiply(vector_to_2d_matrix(next_line[j - 1]),
                            matrix_multiply(binary_exponentiation(T, gap), binary_exponentiation(D, len)));
        // hit = P[j-1, k-1] * T^gap * (T^len - D^l)
        std::vector<std::vector<long double>> hit = matrix_multiply(
            vector_to_2d_matrix(prev_line[j - 1]),
            matrix_multiply(binary_exponentiation(T, gap),
                            subtract_matrices(binary_exponentiation(T, len), binary_exponentiation(D, len))));

        // P[j,k] = dont_hit + hit
        next_line[j] = matrix_to_vector(add_matrices(dont_hit, hit));
      }

      last_col[k][0] = next_line[next_line.size() - 1][0];
      last_col[k][1] = next_line[next_line.size() - 1][1];

      prev_line = next_line;
    }

    std::vector<std::vector<long double>> cur_probs(2, std::vector<long double>(m + 1));
    for (int k = 0; k <= m; k++) {
      // length of gap from end of last interval to end of window
      long long trailing_gap = window_end - ref_intervals_augmented[m].end;
      std::vector<long double> actual_last_col =
          matrix_to_vector(matrix_multiply(vector_to_2d_matrix(last_col[k]), binary_exponentiation(T, trailing_gap)));
      std::cout << "last_col_fast: " << actual_last_col[0] << " " << actual_last_col[1] << "\n";
      for (int ending_state : {0, 1})
        cur_probs[ending_state][k] = log(actual_last_col[ending_state]);
    }
    probs[start_state] = cur_probs;
  }

  std::cout << probs[0][0].size() << " probs.size()\n";

  return probs;
}

std::vector<long double> Model::eval_probs_single_chr_direct_eigen(std::vector<Interval> ref_intervals,
                                                                   std::vector<Interval> query_intervals,
                                                                   const MarkovChain &markov_chain,
                                                                   long long chr_size) {
  logger.error("This is not implemented.");
  exit(1);
}
