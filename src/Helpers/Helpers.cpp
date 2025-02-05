#include "Helpers.h"
#include "../Interval/Interval.h"
#include "../Logger/Logger.h"
#include <algorithm>
#include <fstream>
#include <sstream>
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

std::unordered_map<std::string, long long>
load_chr_sizes(const std::string &file_path) {
  std::unordered_map<std::string, long long> chr_sizes;

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
