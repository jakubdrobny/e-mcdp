#include "WindowModel.hpp"
#include "../WindowResult/WindowResult.hpp"
#include <algorithm>
#include <iostream>

WindowModel::WindowModel() {}
WindowModel::WindowModel(std::vector<Interval> windows,
                         std::vector<Interval> ref_intervals,
                         std::vector<Interval> query_intervals,
                         ChrSizesMap chr_sizes_map, Algorithm algorithm)
    : windows(windows), ref_intervals(ref_intervals),
      query_intervals(query_intervals), algorithm(algorithm) {

  chr_sizes = chr_sizes_map_to_array(chr_sizes_map);
  std::sort(chr_sizes.begin(), chr_sizes.end());
}

// requires the intervals be non-overlapping and windows have non-decreasing
// starts and ends
std::vector<std::vector<Interval>>
WindowModel::get_windows_intervals(const std::vector<Interval> &windows,
                                   const std::vector<Interval> &intervals) {
  std::vector<std::vector<Interval>> results(windows.size());
  if (intervals.empty()) {
    return results;
  }

  std::vector<Interval> sortedIntervals = intervals;
  std::sort(sortedIntervals.begin(), sortedIntervals.end());

  for (size_t idx = 1; idx < sortedIntervals.size(); idx++) {
    if (sortedIntervals[idx].begin < sortedIntervals[idx - 1].end) {
      logger.error("Intervals split into windows need to be non-overlapping.");
      exit(1);
    }
  }

  std::vector<Interval> sortedWindows = windows;
  std::sort(sortedWindows.begin(), sortedWindows.end());

  for (size_t idx = 1; idx < sortedWindows.size(); idx++) {
    if (sortedWindows[idx].end < sortedWindows[idx - 1].end) {
      logger.error("Windows need to have non-decreasing coordinates.");
      exit(1);
    }
  }

  for (size_t windows_idx = 0, intervals_idx = 0;
       windows_idx < sortedWindows.size(); windows_idx++) {
    while (intervals_idx < sortedIntervals.size() &&
           sortedIntervals[intervals_idx].end <= windows[windows_idx].begin)
      intervals_idx++;
    size_t r = intervals_idx;

    std::vector<Interval> cur_result;
    while (r < sortedIntervals.size() &&
           intervals[r].begin < windows[windows_idx].end) {
      Interval sliced_interval = sortedIntervals[r++];
      sliced_interval.begin =
          std::max(sliced_interval.begin, windows[windows_idx].begin);
      sliced_interval.end =
          std::min(sliced_interval.end, windows[windows_idx].end);
      if (sliced_interval.end - sliced_interval.begin < 1)
        continue;
      cur_result.push_back(sliced_interval);
    }
    results[windows_idx] = cur_result;
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

  std::vector<WindowResult> probs_by_window;

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
          count_overlaps_single_chr(ref_intervals_by_window[window_idx],
                                    query_intervals_by_window[window_idx]);
      std::vector<long double> probs = eval_probs_single_chr_direct(
          ref_intervals_by_window[window_idx],
          query_intervals_by_window[window_idx], chr_size);
      Interval cur_window = windows_by_chr[chr_sizes_idx][window_idx];
      probs_by_window.push_back(WindowResult(cur_window, overlap_count, probs));
    }
  }

  return probs_by_window;
}
