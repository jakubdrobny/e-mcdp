#include "WindowModel.hpp"
#include "../Helpers/Helpers.hpp"
#include "../Interval/Section.hpp"
#include "../Results/WindowResult.hpp"
#include <algorithm>
#include <set>

WindowModel::WindowModel() {}
WindowModel::WindowModel(std::vector<Interval> windows, std::vector<Interval> ref_intervals,
                         std::vector<Interval> query_intervals, ChrSizesMap chr_sizes_map, Algorithm algorithm)
    : windows(windows), ref_intervals(ref_intervals), query_intervals(query_intervals), algorithm(algorithm) {

  chr_sizes = chr_sizes_map_to_array(chr_sizes_map);
  std::sort(chr_sizes.begin(), chr_sizes.end());
}

// requires the intervals be non-overlapping
std::vector<std::vector<Interval>> WindowModel::get_windows_intervals_naive(const std::vector<Interval> &windows,
                                                                            const std::vector<Interval> &intervals) {
  std::vector<std::vector<Interval>> results(windows.size());
  if (intervals.empty()) {
    return results;
  }

  std::vector<Interval> sortedIntervals = intervals;
  for (size_t idx = 0; idx + 1 < sortedIntervals.size(); idx++) {
    if (sortedIntervals[idx].end > sortedIntervals[idx + 1].begin) {
      logger.error("intervals must be non-overlapping for splitting into windows");
      exit(1);
    }
  }

  for (size_t windows_idx = 0; windows_idx < windows.size(); windows_idx++) {
    Interval window = windows[windows_idx];
    for (Interval interval : intervals) {
      if (interval.end <= window.begin || interval.begin >= window.end)
        continue;
      Interval sliced_interval = interval;
      sliced_interval.begin = std::max(sliced_interval.begin, window.begin);
      sliced_interval.end = std::min(sliced_interval.end, window.end);
      if (sliced_interval.end - sliced_interval.begin < 1)
        continue;
      results[windows_idx].push_back(sliced_interval);
    }
  }

  return results;
}

template <typename WindowType>
std::vector<std::vector<Interval>> WindowModel::get_windows_intervals(const std::vector<WindowType> &windows,
                                                                      const std::vector<Interval> &intervals) {
  static_assert(std::is_base_of<Interval, WindowType>::value, "WindowType must inherit from Interval");

  std::vector<std::vector<Interval>> results(windows.size());
  if (intervals.empty()) {
    return results;
  }

  if (!are_intervals_non_overlapping(intervals)) {
    logger.error("intervals need to be non-overlapping for spliting into "
                 "windows to happen.");
    exit(1);
  }

  struct Event {
    long long pos;
    bool start, interval;
    size_t idx;
  };

  std::vector<Event> events((windows.size() << 1) + (intervals.size() << 1));

  for (size_t idx = 0; idx < windows.size(); idx++) {
    Interval window = windows[idx];
    events[idx << 1] = {window.get_begin(), 1, 0, idx};
    events[(idx << 1) | 1] = {window.get_end(), 0, 0, idx};
  }

  size_t buf = windows.size() << 1;
  for (size_t idx = 0; idx < intervals.size(); idx++) {
    Interval interval = intervals[idx];
    events[buf + (idx << 1)] = {interval.begin, 1, 1, idx};
    events[buf + (idx << 1) + 1] = {interval.end, 0, 1, idx};
  }

  std::sort(events.begin(), events.end(), [](const Event &e1, const Event &e2) {
    if (e1.pos == e2.pos) {
      if (e1.start == e2.start) {
        if (e1.interval == e2.interval) {
          return e1.idx < e2.idx;
        }
        return e1.interval < e2.interval;
      }
      return e1.start < e2.start;
    }
    return e1.pos < e2.pos;
  });

  std::set<int> opened_intervals, opened_windows;
  for (Event event : events) {
    if (event.start) {
      if (event.interval) {
        opened_intervals.insert(event.idx);
      } else {
        opened_windows.insert(event.idx);
      }
    } else {
      // ends of windows should come first so that we
      // do not close the intervals before closing the windows
      if (!event.interval) {
        for (int interval_idx : opened_intervals) {
          Interval sliced_interval = slice_interval_by_window(windows[event.idx], intervals[interval_idx]);
          if (sliced_interval.length() < 1)
            continue;
          results[event.idx].push_back(sliced_interval);
        }
        opened_windows.erase(event.idx);
      } else {
        // if there are windows opened add this interval to them
        for (int window_idx : opened_windows) {
          Interval sliced_interval = slice_interval_by_window(windows[window_idx], intervals[event.idx]);
          if (sliced_interval.length() < 1)
            continue;
          results[window_idx].push_back(sliced_interval);
        }
        opened_intervals.erase(event.idx);
      }
    }
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

  std::vector<std::vector<Interval>> windows_by_chr(chr_sizes.size()), ref_intervals_by_chr(chr_sizes.size()),
      query_intervals_by_chr(chr_sizes.size());

  for (size_t chr_sizes_idx = 0, windows_idx = 0, ref_idx = 0, query_idx = 0; chr_sizes_idx < chr_sizes.size();
       chr_sizes_idx++) {
    std::string chr_name = chr_sizes[chr_sizes_idx].first;

    windows_by_chr[chr_sizes_idx] = select_intervals_by_chr_name(windows, windows_idx, chr_name);
    ref_intervals_by_chr[chr_sizes_idx] = select_intervals_by_chr_name(ref_intervals, ref_idx, chr_name),
    query_intervals_by_chr[chr_sizes_idx] = select_intervals_by_chr_name(query_intervals, query_idx, chr_name);
  }

  std::vector<WindowResult> probs_by_window;

  // turn off for debugging
  // #pragma omp parallel for
  for (size_t chr_sizes_idx = 0; chr_sizes_idx < chr_sizes.size(); chr_sizes_idx++) {
    std::vector<WindowResult> chromosome_probs_by_window;
    if (algorithm == Algorithm::NAIVE) {
      chromosome_probs_by_window =
          probs_by_window_single_chr_naive(windows_by_chr[chr_sizes_idx], ref_intervals_by_chr[chr_sizes_idx],
                                           query_intervals_by_chr[chr_sizes_idx], chr_sizes[chr_sizes_idx]);
    } else if (algorithm == Algorithm::FAST) {
      chromosome_probs_by_window =
          probs_by_window_single_chr_smarter(windows_by_chr[chr_sizes_idx], ref_intervals_by_chr[chr_sizes_idx],
                                             query_intervals_by_chr[chr_sizes_idx], chr_sizes[chr_sizes_idx]);
    } else {
      logger.error("invalid algorithm.");
      exit(1);
    }

    probs_by_window.reserve(probs_by_window.size() + chromosome_probs_by_window.size());
    probs_by_window.insert(probs_by_window.end(), chromosome_probs_by_window.begin(), chromosome_probs_by_window.end());
  }

  return probs_by_window;
}

std::vector<WindowResult> WindowModel::probs_by_window_single_chr_naive(
    const std::vector<Interval> &windows, const std::vector<Interval> &ref_intervals,
    const std::vector<Interval> &query_intervals, const std::pair<std::string, long long> chr_size_entry) {

  std::string chr_name = chr_size_entry.first;
  logger.info("Loading windows and their intervals for chromosome: " + chr_name);

  long long chr_size = chr_size_entry.second;

  std::vector<std::vector<Interval>> ref_intervals_by_window = get_windows_intervals<Interval>(windows, ref_intervals);
  std::vector<std::vector<Interval>> query_intervals_by_window =
      get_windows_intervals<Interval>(windows, query_intervals);

  logger.info("Calculating probs for windows in chromsome: " + chr_name);

  std::vector<WindowResult> probs_by_window;

  MarkovChain markov_chain(chr_size, query_intervals);

  // TODO: can this work?
  // #pragma omp parallel for
  for (size_t window_idx = 0; window_idx < windows.size(); window_idx++) {
    long long overlap_count =
        count_overlaps_single_chr(ref_intervals_by_window[window_idx], query_intervals_by_window[window_idx]);
    std::vector<long double> probs = eval_probs_single_chr_direct(
        ref_intervals_by_window[window_idx], query_intervals_by_window[window_idx], markov_chain, chr_size);
    Interval cur_window = windows[window_idx];
    probs_by_window.push_back(WindowResult(cur_window, overlap_count, probs));
  }

  return probs_by_window;
}

std::vector<WindowResult> WindowModel::probs_by_window_single_chr_smarter(
    const std::vector<Interval> &windows, const std::vector<Interval> &ref_intervals,
    const std::vector<Interval> &query_intervals, const std::pair<std::string, long long> chr_size_entry) {
  if (windows.empty()) {
    return {};
  }

  long long chr_size = chr_size_entry.second;

  // 1. create sections from (possibly) overlapping set of windows
  WindowSectionSplitResult windowSectionSplitResult =
      split_windows_into_non_overlapping_sections(windows, ref_intervals);
  std::vector<Section> sections = windowSectionSplitResult.get_sections();
  std::vector<Interval> spans = windowSectionSplitResult.get_spans();

  // 2. load intervals into sections, will be fast since both are
  // non-overlapping
  std::vector<std::vector<Interval>> ref_intervals_by_section = get_windows_intervals<Section>(sections, ref_intervals),
                                     query_intervals_by_section =
                                         get_windows_intervals<Section>(sections, query_intervals);

  // 3. calculate transition matrices
  MarkovChain markov_chain(chr_size, query_intervals);
  // markov_chain.print();

  // 4. calculature probs of each section
  for (size_t section_idx = 0; section_idx < sections.size(); section_idx++) {
    SectionProbs probs = eval_probs_single_section(ref_intervals_by_section[section_idx], sections[section_idx].begin,
                                                   sections[section_idx].end, markov_chain);
    sections[section_idx].set_probs(probs);
  }

  // 5. merge section probs for each window
  std::vector<WindowResult> probs_by_window(windows.size());

  for (size_t windows_idx = 0; windows_idx < windows.size(); windows_idx++) {
    Interval span = spans[windows_idx];
    Section section = sections[span.begin];

    // merge probs for sections
    for (long long sections_idx = span.begin + 1; sections_idx < span.end; sections_idx++) {
      section = join_sections(section, sections[sections_idx], markov_chain);
    }

    // merge the final 4 sets of probs for window into one
    std::vector<long double> cur_windows_single_probs =
        merge_multi_probs(section.get_probs().get_normal(), markov_chain);
    probs_by_window[windows_idx] = WindowResult(windows[windows_idx], 01, cur_windows_single_probs);
  }

  return probs_by_window;
}

SectionProbs WindowModel::eval_probs_single_section(const std::vector<Interval> &ref_intervals, long long section_start,
                                                    long long section_end, const MarkovChain &markov_chain) {
  MultiProbs probs_normal = eval_probs_single_chr_direct_new(ref_intervals, section_start, section_end, markov_chain);

  std::vector<Interval> ref_intervals_except_last(ref_intervals.begin(),
                                                  ref_intervals.end() - (!ref_intervals.empty()));
  long long new_section_end = ref_intervals_except_last.empty() ? section_start : ref_intervals_except_last.back().end;
  MultiProbs probs_except_last =
      eval_probs_single_chr_direct_new(ref_intervals_except_last, section_start, new_section_end, markov_chain);

  std::vector<Interval> ref_intervals_except_first(ref_intervals.begin() + (!ref_intervals.empty()),
                                                   ref_intervals.end());
  long long new_section_start = ref_intervals.empty() ? section_start : ref_intervals[0].end;
  MultiProbs probs_except_first =
      eval_probs_single_chr_direct_new(ref_intervals_except_first, new_section_start, section_end, markov_chain);

  std::vector<Interval> ref_intervals_except_first_and_last(
      ref_intervals_except_first.begin(), ref_intervals_except_first.end() - (!ref_intervals_except_first.empty()));
  MultiProbs probs_except_first_and_last = eval_probs_single_chr_direct_new(
      ref_intervals_except_first_and_last, new_section_start, new_section_end, markov_chain);

  return SectionProbs(probs_normal, probs_except_first, probs_except_last, probs_except_first_and_last);
}

template std::vector<std::vector<Interval>>
WindowModel::get_windows_intervals<Interval>(const std::vector<Interval> &windows,
                                             const std::vector<Interval> &intervals);

template std::vector<std::vector<Interval>>
WindowModel::get_windows_intervals<Section>(const std::vector<Section> &sections,
                                            const std::vector<Interval> &intervals);
