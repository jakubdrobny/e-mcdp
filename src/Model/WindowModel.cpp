#include "WindowModel.hpp"
#include "../WindowResult/WindowResult.hpp"
#include <algorithm>

WindowModel::WindowModel() {}
WindowModel::WindowModel(std::vector<Interval> windows,
                         std::vector<Interval> ref_intervals,
                         std::vector<Interval> query_intervals,
                         ChrSizesMap chr_sizes_map)
    : windows(windows), ref_intervals(ref_intervals),
      query_intervals(query_intervals) {

  chr_sizes = chr_sizes_map_to_array(chr_sizes_map);
  std::sort(chr_sizes.begin(), chr_sizes.end());
}

std::vector<std::vector<Interval>>
WindowModel::get_windows_intervals(const std::vector<Interval> &windows,
                                   const std::vector<Interval> &intervals) {
  std::vector<std::vector<Interval>> results(windows.size());
  if (intervals.empty()) {
    return results;
  }

  std::vector<Interval> sortedIntervals = intervals;
  std::sort(sortedIntervals.begin(), sortedIntervals.end());

  std::vector<std::pair<long long, int>> B;
  for (size_t i = 0; i < sortedIntervals.size(); i++) {
    B.emplace_back(sortedIntervals[i].end, i);
  }

  int n = B.size();
  std::vector<std::vector<std::pair<long long, int>>> BIT(n + 1);

  for (int i = 0; i < n; i++) {
    int idx = i + 1;
    while (idx <= n) {
      BIT[idx].push_back(B[i]);
      idx += idx & -idx;
    }
  }

  for (auto &list : BIT) {
    sort(list.begin(), list.end());
  }

  const std::string chr_name = sortedIntervals[0].chr_name;
  for (size_t windows_idx = 0; windows_idx < windows.size(); windows_idx++) {
    const Interval window = windows[windows_idx];
    int k = std::upper_bound(sortedIntervals.begin(), sortedIntervals.end(),
                             Interval(chr_name, window.end, 0),
                             [](const Interval &a, const Interval &b) {
                               return a.begin < b.begin;
                             }) -
            sortedIntervals.begin();

    std::vector<int> indices;
    int idx = k;
    while (idx > 0) {
      auto &list = BIT[idx];
      auto it = std::upper_bound(list.begin(), list.end(),
                                 std::make_pair(window.begin, 0));
      for (auto iter = it; iter != list.end(); iter++) {
        indices.push_back(iter->second);
      }
      idx -= idx & -idx;
    }

    sort(indices.begin(), indices.end());
    indices.erase(unique(indices.begin(), indices.end()), indices.end());

    std::vector<Interval> overlapping;
    for (int i : indices) {
      overlapping.push_back(
          Interval(chr_name, std::max(window.begin, sortedIntervals[i].begin),
                   std::min(window.end, sortedIntervals[i].end)));
    }

    results[windows_idx] = overlapping;
  }

  return results;
}

std::vector<WindowResult> WindowModel::run() {
  logger.info("Running WindowModel...");
  logger.info("Sorting intervals and windows...");

  std::sort(ref_intervals.begin(), ref_intervals.end());
  std::sort(query_intervals.begin(), query_intervals.end());
  std::sort(windows.begin(), windows.end());

  logger.info("Grouping intervals and windows by chromosome...");

  std::vector<std::vector<Interval>> windows_by_chr(chr_sizes.size()),
      ref_intervals_by_chr(chr_sizes.size()),
      query_intervals_by_chr(chr_sizes.size());

  for (size_t chr_sizes_idx = 0, windows_idx = 0, ref_idx = 0, query_idx = 0;
       chr_sizes_idx < chr_sizes.size(); chr_sizes_idx++) {
    std::string chr_name = chr_sizes[chr_sizes_idx].first;

    windows_by_chr[chr_sizes_idx] =
        select_intervals_by_chr_name(windows, windows_idx, chr_name);
    ref_intervals_by_chr[chr_sizes_idx] =
        select_intervals_by_chr_name(ref_intervals, ref_idx, chr_name),
    query_intervals_by_chr[chr_sizes_idx] =
        select_intervals_by_chr_name(query_intervals, query_idx, chr_name);
  }

  std::vector<WindowResult> probs_by_window(windows.size());

  // turn off for debugging
  // #pragma omp parallel for
  for (size_t chr_sizes_idx = 0; chr_sizes_idx < chr_sizes.size();
       chr_sizes_idx++) {
    std::string chr_name = chr_sizes[chr_sizes_idx].first;
    logger.info("Loading windows and their intervals for chromosome: " +
                chr_name);

    long long chr_size = chr_sizes[chr_sizes_idx].second;

    std::vector<std::vector<Interval>> ref_intervals_by_window =
        get_windows_intervals(windows_by_chr[chr_sizes_idx],
                              ref_intervals_by_chr[chr_sizes_idx]);
    std::vector<std::vector<Interval>> query_intervals_by_window =
        get_windows_intervals(windows_by_chr[chr_sizes_idx],
                              query_intervals_by_chr[chr_sizes_idx]);

    logger.info("Calculating probs for windows in chromsome: " + chr_name);
    // TODO: can this work?
    // #pragma omp parallel for
    for (size_t window_idx = 0;
         window_idx < windows_by_chr[chr_sizes_idx].size(); window_idx++) {
      long long overlap_count =
          count_overlaps(ref_intervals_by_window[window_idx],
                         query_intervals_by_window[window_idx]);
      std::vector<long double> probs = eval_probs_single_chr_direct(
          ref_intervals_by_chr[chr_sizes_idx],
          query_intervals_by_chr[chr_sizes_idx], chr_size);
      probs_by_window.push_back(WindowResult(
          windows_by_chr[chr_sizes_idx][window_idx], overlap_count, probs));
    }
  }

  return probs_by_window;
}
