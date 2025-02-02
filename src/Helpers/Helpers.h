#ifndef HELPERS_H
#define HELPERS_H

#include "../Interval/Interval.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

Interval parse_intervals_line(std::string line);
std::vector<Interval> load_intervals(const std::string &file_path,
                                     bool is_closed = false);
std::unordered_map<std::string, long long>
load_chr_sizes(const std::string &file_path);
std::unordered_set<std::string> load_chr_names_from_chr_sizes(
    const std::unordered_map<std::string, long long int> &chr_sizes);
std::vector<Interval>
filter_intervals_by_chr_name(std::vector<Interval> intervals,
                             std::unordered_set<std::string> chr_names);
std::vector<Interval>
merge_non_disjoint_intervals(std::vector<Interval> intervals);
std::vector<Interval> remove_empty_intervals(std::vector<Interval> intervals);

#endif
