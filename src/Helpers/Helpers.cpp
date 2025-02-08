#include "Helpers.h"
#include "../Interval/Interval.h"
#include "../Logger/Logger.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// if exp_elem_cnt = 0, then parses any number of values
// otherwise if parsed values != exp_elem_cnt, exits
std::vector<std::string> split_string(std::string str, char delimeter,
                                      size_t exp_elem_cnt) {
  std::vector<std::string> vals;

  std::stringstream ss(str);
  std::string word;
  while (getline(ss, word, delimeter)) {
    vals.push_back(word);
    if (exp_elem_cnt > 0 && vals.size() > exp_elem_cnt)
      break;
  }

  if (exp_elem_cnt > 0 && vals.size() != exp_elem_cnt) {
    logger.error("Expected " + std::to_string(exp_elem_cnt) + " values in " +
                 str + " split by " + delimeter + ". Exiting...");
    exit(1);
  }

  return vals;
}

// expects {chr_name} {begin} {end}
Interval parse_intervals_line(std::string line) {
  size_t exp_elem_cnt = 3;
  std::vector<std::string> vals = split_string(line, '\t', exp_elem_cnt);

  if (vals.size() != exp_elem_cnt) {
    logger.error("Invalid line format on line: " + line +
                 ". Should be {chr_name} {begin} {end}. Exiting...");
    exit(1);
  }

  return Interval(vals[0], std::stoll(vals[1]), std::stoll(vals[2]));
}

std::vector<Interval> load_intervals(const std::string &file_path,
                                     bool is_closed) {
  std::vector<Interval> intervals;

  std::ifstream input_file(file_path);
  std::string line;
  while (std::getline(input_file, line)) {
    Interval interval = parse_intervals_line(line);
    if (is_closed)
      interval.end--;

    if (interval.begin >= interval.end) {
      logger.error("Begin should be strictly smaller than end in stated "
                   "intervals. Exiting...");
      exit(1);
    }

    if (interval.begin < 0 || interval.end < 0) {
      logger.error("Interval bounds should be non-negative. Exiting...");
      exit(1);
    }

    intervals.push_back(interval);
  }

  return intervals;
}

ChrSizesMap load_chr_sizes(const std::string &file_path) {
  ChrSizesMap chr_sizes;

  std::ifstream input_file(file_path);
  std::string line;
  while (getline(input_file, line)) {
    std::vector<std::string> vals = split_string(line, '\t', 2);
    std::string chr_name = vals[0];
    long long chr_size = std::stoll(vals[1]);

    if (chr_sizes.count(chr_name)) {
      logger.error("Size for chromosome " + chr_name +
                   " specified more than once. Exiting...");
      exit(1);
    }

    if (chr_size < 0) {
      logger.error("Chromosome sizes should be position numbers. Exiting...");
      exit(1);
    }

    chr_sizes[chr_name] = chr_size;
  }

  return chr_sizes;
}

std::unordered_set<std::string> load_chr_names_from_chr_sizes(
    const std::unordered_map<std::string, long long> &chr_sizes) {
  std::unordered_set<std::string> chr_names;
  for (const auto &p : chr_sizes) {
    chr_names.insert(p.first);
  }
  return chr_names;
}

std::vector<Interval>
filter_intervals_by_chr_name(std::vector<Interval> intervals,
                             std::unordered_set<std::string> chr_names) {
  std::vector<Interval> new_intervals;

  for (Interval interval : intervals) {
    if (chr_names.count(interval.chr_name)) {
      new_intervals.push_back(interval);
    }
  }

  return new_intervals;
}

std::vector<Interval>
merge_non_disjoint_intervals(std::vector<Interval> intervals) {
  if (intervals.size() == 0) {
    return intervals;
  }

  std::sort(intervals.begin(), intervals.end());

  std::vector<Interval> new_intervals;
  Interval cur_interval = intervals[0];
  for (size_t interval_idx = 1; interval_idx < intervals.size();
       interval_idx++) {
    Interval new_interval = intervals[interval_idx];
    if (cur_interval.chr_name != new_interval.chr_name ||
        cur_interval.end < new_interval.begin) {
      new_intervals.push_back(cur_interval);
      cur_interval = new_interval;
      continue;
    }

    if (new_interval.end > cur_interval.end) {
      cur_interval.end = new_interval.end;
    }
  }

  new_intervals.push_back(cur_interval);

  return new_intervals;
}

std::vector<Interval> remove_empty_intervals(std::vector<Interval> intervals) {
  std::vector<Interval> new_intervals;

  for (Interval interval : intervals) {
    if (interval.end - interval.begin > 0) {
      new_intervals.push_back(interval);
    }
  }

  return new_intervals;
}

std::vector<std::string>
get_sorted_chr_names_from_intervals(std::vector<Interval> intervals) {
  std::vector<std::string> chr_names;

  if (intervals.empty()) {
    return chr_names;
  }

  sort(intervals.begin(), intervals.end());
  chr_names.push_back(intervals[0].chr_name);
  for (Interval interval : intervals) {
    if (interval.chr_name != intervals.back().chr_name) {
      chr_names.push_back(interval.chr_name);
    }
  }

  return chr_names;
}

template <typename T>
void extend(std::vector<T> &self, const std::vector<T> &other) {
  self.reserve(self.size() + other.size());
  for (const auto &element : other) {
    self.push_back(element);
  }
}

// assuming intervals are sorted
long long count_overlaps_single_chr(std::vector<Interval> ref_intervals,
                                    std::vector<Interval> query_intervals) {
  long long overlap_count = 0;
  bool is_ref_interval_open = false, is_query_interval_open = false,
       is_current_ref_interval_counted = false;

  std::vector<std::vector<long long>> events;
  for (Interval interval : ref_intervals) {
    events.push_back({interval.begin, 0, 0});
    events.push_back({interval.end, 0, 1});
  }
  for (Interval interval : query_intervals) {
    events.push_back({interval.begin, 1, 0});
    events.push_back({interval.end, 1, 1});
  }
  sort(events.begin(), events.end());

  long long last_pos = -1;
  for (std::vector<long long> event : events) {
    long long pos = event[0];
    bool is_query = event[1], is_end = event[2];
    if (last_pos < pos) {
      last_pos = pos;
      if (is_ref_interval_open && is_query_interval_open &&
          !is_current_ref_interval_counted) {
        overlap_count++;
        is_current_ref_interval_counted = true;
      }
    }
    if (!is_query && !is_end) {
      is_ref_interval_open = true;
      is_current_ref_interval_counted = false;
    }
    if (!is_query && is_end) {
      is_ref_interval_open = false;
    }
    if (is_query && !is_end) {
      is_query_interval_open = true;
    }
    if (is_query && is_end) {
      is_query_interval_open = false;
    }
  }

  return overlap_count;
}

long long count_overlaps(std::vector<Interval> ref_intervals,
                         std::vector<Interval> query_intervals) {
  std::sort(ref_intervals.begin(), ref_intervals.end());
  std::sort(query_intervals.begin(), query_intervals.end());

  std::vector<std::string> chr_names,
      ref_chr_names = get_sorted_chr_names_from_intervals(ref_intervals),
      query_chr_names = get_sorted_chr_names_from_intervals(query_intervals);
  extend(chr_names, ref_chr_names);
  extend(chr_names, query_chr_names);

  auto get_chr_intervals = [](std::vector<Interval> &intervals, int &idx,
                              std::string chr_name) {
    std::vector<Interval> chr_intervals;
    for (; intervals[idx].chr_name <= chr_name; idx++)
      if (intervals[idx].chr_name == chr_name)
        chr_intervals.push_back(intervals[idx]);
    return chr_intervals;
  };

  int ref_idx = 0, query_idx = 0;
  long long total_overlap_count = 0;
  for (std::string chr_name : chr_names) {
    std::vector<Interval> chr_ref_intervals = get_chr_intervals(
                              ref_intervals, ref_idx, chr_name),
                          chr_query_intervals = get_chr_intervals(
                              query_intervals, query_idx, chr_name);

    if (chr_ref_intervals.empty() || chr_query_intervals.empty())
      continue;

    long long chr_overlap_count =
        count_overlaps_single_chr(chr_ref_intervals, chr_query_intervals);
    total_overlap_count += chr_overlap_count;
  }

  return total_overlap_count;
}

ChrSizesVector chr_sizes_map_to_array(ChrSizesMap &chr_sizes_map) {
  ChrSizesVector chr_sizes_vector;
  for (std::pair<std::string, long long> p : chr_sizes_map)
    chr_sizes_vector.push_back(p);
  sort(chr_sizes_vector.begin(), chr_sizes_vector.end());
  return chr_sizes_vector;
}

long double calculate_joint_pvalue(
    std::vector<std::vector<long double>> &probs_by_chromosome,
    long long overlap_count) {
  return 0.;
}
