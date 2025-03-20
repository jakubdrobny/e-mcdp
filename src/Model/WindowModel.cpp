#include "WindowModel.hpp"
#include "../WindowResult/WindowResult.hpp"
#include <algorithm>
#include <iostream>
#include <set>

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

// requires the intervals be non-overlapping
std::vector<std::vector<Interval>>
WindowModel::get_windows_intervals(const std::vector<Interval> &windows,
                                   const std::vector<Interval> &intervals) {
  std::vector<std::vector<Interval>> results(windows.size());
  if (intervals.empty()) {
    return results;
  }

  // sorted by end increasing primarily and begin secondarily
  std::vector<Interval> sortedIntervals = intervals;
  std::sort(sortedIntervals.begin(), sortedIntervals.end());

  for (size_t idx = 1; idx < sortedIntervals.size(); idx++) {
    if (sortedIntervals[idx].begin < sortedIntervals[idx - 1].end) {
      logger.error("Intervals split into windows need to be non-overlapping.");
      exit(1);
    }
  }

  using WindowEntry = std::pair<Interval, int>;
  std::vector<WindowEntry> sortedWindowsByBegin(windows.size());
  std::vector<WindowEntry> sortedWindowsByEnd(windows.size());
  for (size_t windows_idx = 0; windows_idx < windows.size(); windows_idx++) {
    sortedWindowsByBegin[windows_idx] = {windows[windows_idx], windows_idx};
    sortedWindowsByEnd[windows_idx] = {windows[windows_idx], windows_idx};
  }

  std::sort(sortedWindowsByBegin.begin(), sortedWindowsByBegin.end());
  std::sort(sortedWindowsByEnd.begin(), sortedWindowsByEnd.end(),
            [](const WindowEntry &e1, const WindowEntry &e2) {
              Interval i1 = e1.first, i2 = e2.first;
              if (i1.end == i2.end) {
                if (i1.begin == i2.begin)
                  return i1.chr_name < i2.chr_name;
                return i1.begin < i2.begin;
              }
              return i1.end < i2.end;
            });

  std::vector<std::set<size_t>> resultsSet(windows.size());
  for (size_t intervals_idx = 0; intervals_idx < sortedIntervals.size();
       intervals_idx++) {
    Interval interval = sortedIntervals[intervals_idx];
    std::cout << "interval: " << interval << "\n";
    auto find_intervals = [interval, intervals_idx, &resultsSet](
                              const std::vector<WindowEntry> &sortedWindows) {
      size_t lb =
          std::lower_bound(sortedWindows.begin(), sortedWindows.end(), interval,
                           [](const WindowEntry &e1, const Interval &i2) {
                             return e1.first.begin < i2.begin;
                           }) -
          sortedWindows.begin();
      size_t ub =
          std::lower_bound(sortedWindows.begin(), sortedWindows.end(), interval,
                           [](const WindowEntry &e1, const Interval &i2) {
                             return e1.first.begin < i2.end;
                           }) -
          sortedWindows.begin();
      std::cout << "lb: " << lb << ", ub: " << ub << ", sortedWindows:";
      for (WindowEntry _we : sortedWindows) {
        std::cout << " " << _we.first;
      }
      std::cout << "\n";

      for (size_t idx = lb; idx < ub; idx++) {
        WindowEntry we = sortedWindows[idx];
        Interval sliced_interval = interval, window = we.first;
        sliced_interval.begin = std::max(sliced_interval.begin, window.begin);
        sliced_interval.end = std::min(sliced_interval.end, window.end);
        if (sliced_interval.end - sliced_interval.begin < 1)
          continue;
        resultsSet[we.second].insert(intervals_idx);
      }
    };

    find_intervals(sortedWindowsByBegin);
    find_intervals(sortedWindowsByEnd);
  }

  for (size_t idx = 0; idx < resultsSet.size(); idx++) {
    std::cout << "window: " << windows[idx] << ":";
    for (size_t intervals_idx : resultsSet[idx]) {
      std::cout << " " << sortedIntervals[intervals_idx];
      results[idx].push_back(sortedIntervals[intervals_idx]);
    }
    std::cout << "\n";
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
