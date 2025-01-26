#ifndef HELPERS_H 
#define HELPERS_H

#include "../Interval/Interval.h"

#include <string>
#include <unordered_map>
#include <vector>

Interval parse_intervals_line(std::string line);
std::vector<Interval> load_intervals(const std::string& file_path, bool is_closed = false);
std::unordered_map<std::string, long long > load_chr_sizes(const std::string& file_path);

#endif
