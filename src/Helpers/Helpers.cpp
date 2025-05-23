#include "Helpers.hpp"
#include "../Interval/Interval.hpp"
#include "../Interval/Section.hpp"
#include "../Logger/Logger.hpp"
#include "../Model/Model.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// if exp_elem_cnt = 0, then parses any number of values
// otherwise if parsed values != exp_elem_cnt, exits
std::vector<std::string> split_string(std::string str, char delimeter, size_t exp_elem_cnt) {
  std::vector<std::string> vals;

  std::stringstream ss(str);
  std::string word;
  while (getline(ss, word, delimeter)) {
    vals.push_back(word);
    if (exp_elem_cnt > 0 && vals.size() > exp_elem_cnt)
      break;
  }

  if (exp_elem_cnt > 0 && vals.size() != exp_elem_cnt) {
    logger.error("Expected " + std::to_string(exp_elem_cnt) + " values in " + str + " split by " + delimeter + ".");
    exit(1);
  }

  return vals;
}

// expects {chr_name} {begin} {end}
Interval parse_intervals_line(std::string line) {
  size_t exp_elem_cnt = 3;
  std::vector<std::string> vals = split_string(line, '\t', exp_elem_cnt);

  if (vals.size() != exp_elem_cnt) {
    logger.error("Invalid line format on line: " + line + ". Should be {chr_name} {begin} {end}.");
    exit(1);
  }

  return Interval(vals[0], std::stoll(vals[1]), std::stoll(vals[2]));
}

std::vector<Interval> load_intervals(const std::string &file_path, bool is_closed) {
  std::vector<Interval> intervals;

  std::ifstream input_file(file_path);
  if (!input_file.is_open()) {
    logger.error("Failed to open intervals file: " + file_path);
    exit(1);
  }

  std::string line;
  while (std::getline(input_file, line)) {
    Interval interval = parse_intervals_line(line);
    if (is_closed)
      interval.end--;

    if (interval.begin >= interval.end) {
      logger.error("Begin should be strictly smaller than end in stated "
                   "intervals.");
      exit(1);
    }

    if (interval.begin < 0 || interval.end < 0) {
      logger.error("Interval bounds should be non-negative.");
      exit(1);
    }

    intervals.push_back(interval);
  }

  return intervals;
}

// assumes args.check_invalid_args has already been run
std::vector<Interval> load_windows(Args &args, std::unordered_map<std::string, long long> &chr_sizes) {
  std::vector<Interval> windows;

  if (args.windows_source == "basic" || args.windows_source == "dense") {
    std::vector<std::pair<std::string, long long>> chr_sizes_vec(chr_sizes.begin(), chr_sizes.end());
    std::sort(chr_sizes_vec.begin(), chr_sizes_vec.end());
    for (std::pair<std::string, long long> chr : chr_sizes) {
      std::string chr_name = chr.first;
      long long chr_size = chr.second;
      long long l = 0, r = std::min(chr_size, args.windows_size);
      while (1) {
        windows.push_back({chr_name, l, r});
        if (r >= chr_size) {
          break;
        }

        if (args.windows_source == "basic") {
          l = r;
          r = std::min(chr_size, r + args.windows_size);
        } else if (args.windows_source == "dense") {
          l += args.windows_step;
          r = std::min(chr_size, r + args.windows_step);
        } else {
          logger.error("If you see this, something went horribly wrong :D");
          exit(1);
        }
      }
    }
  } else if (args.windows_source == "file") {
    windows = load_intervals(args.windows_path);
  }

  return windows;
}

ChrSizesMap load_chr_sizes(const std::string &file_path) {
  ChrSizesMap chr_sizes;

  std::ifstream input_file(file_path);
  if (!input_file.is_open()) {
    logger.error("Failed to open chromosome sizes file: " + file_path);
    exit(1);
  }

  std::string line;
  while (getline(input_file, line)) {
    std::vector<std::string> vals = split_string(line, '\t', 2);
    std::string chr_name = vals[0];
    long long chr_size = std::stoll(vals[1]);

    if (chr_sizes.count(chr_name)) {
      logger.error("Size for chromosome " + chr_name + " specified more than once.");
      exit(1);
    }

    if (chr_size < 0) {
      logger.error("Chromosome sizes should be position numbers.");
      exit(1);
    }

    chr_sizes[chr_name] = chr_size;
  }

  return chr_sizes;
}

std::unordered_set<std::string>
load_chr_names_from_chr_sizes(const std::unordered_map<std::string, long long> &chr_sizes) {
  std::unordered_set<std::string> chr_names;
  for (const auto &p : chr_sizes) {
    chr_names.insert(p.first);
  }
  return chr_names;
}

std::vector<Interval> filter_intervals_by_chr_name(std::vector<Interval> intervals,
                                                   std::unordered_set<std::string> chr_names) {
  std::vector<Interval> new_intervals;

  for (Interval interval : intervals) {
    if (chr_names.count(interval.chr_name)) {
      new_intervals.push_back(interval);
    }
  }

  return new_intervals;
}

std::vector<Interval> merge_non_disjoint_intervals(std::vector<Interval> intervals) {
  if (intervals.size() == 0) {
    return intervals;
  }

  std::sort(intervals.begin(), intervals.end());

  std::vector<Interval> new_intervals;
  Interval cur_interval = intervals[0];
  for (size_t interval_idx = 1; interval_idx < intervals.size(); interval_idx++) {
    Interval new_interval = intervals[interval_idx];
    if (cur_interval.chr_name != new_interval.chr_name || cur_interval.end < new_interval.begin) {
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

std::vector<std::string> get_sorted_chr_names_from_intervals(std::vector<Interval> intervals) {
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

template <typename T> void extend(std::vector<T> &self, const std::vector<T> &other) {
  self.reserve(self.size() + other.size());
  for (const auto &element : other) {
    self.push_back(element);
  }
}

// assuming intervals are sorted
long long count_overlaps_single_chr(std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals) {
  long long overlap_count = 0;
  bool is_ref_interval_open = false, is_query_interval_open = false, is_current_ref_interval_counted = false;

  std::vector<std::vector<long long>> events;
  for (Interval interval : ref_intervals) {
    events.push_back({interval.begin, 0, 1});
    events.push_back({interval.end, 0, 0});
  }
  for (Interval interval : query_intervals) {
    events.push_back({interval.begin, 1, 1});
    events.push_back({interval.end, 1, 0});
  }
  sort(events.begin(), events.end());

  long long last_pos = -1;
  for (std::vector<long long> event : events) {
    long long pos = event[0];
    bool is_query = event[1], is_end = !event[2];
    if (last_pos < pos) {
      last_pos = pos;
      if (is_ref_interval_open && is_query_interval_open && !is_current_ref_interval_counted) {
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

long long count_overlaps(std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals) {
  std::sort(ref_intervals.begin(), ref_intervals.end());
  std::sort(query_intervals.begin(), query_intervals.end());

  std::vector<std::string> chr_names, ref_chr_names = get_sorted_chr_names_from_intervals(ref_intervals),
                                      query_chr_names = get_sorted_chr_names_from_intervals(query_intervals);
  extend(chr_names, ref_chr_names);
  extend(chr_names, query_chr_names);

  auto get_chr_intervals = [](std::vector<Interval> &intervals, int &idx, std::string chr_name) {
    std::vector<Interval> chr_intervals;
    for (; idx < (int)intervals.size() && intervals[idx].chr_name <= chr_name; idx++)
      if (intervals[idx].chr_name == chr_name)
        chr_intervals.push_back(intervals[idx]);
    return chr_intervals;
  };

  int ref_idx = 0, query_idx = 0;
  long long total_overlap_count = 0;
  for (std::string chr_name : chr_names) {
    std::vector<Interval> chr_ref_intervals = get_chr_intervals(ref_intervals, ref_idx, chr_name),
                          chr_query_intervals = get_chr_intervals(query_intervals, query_idx, chr_name);

    if (chr_ref_intervals.empty() || chr_query_intervals.empty())
      continue;

    long long chr_overlap_count = count_overlaps_single_chr(chr_ref_intervals, chr_query_intervals);
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

template void extend<Interval>(std::vector<Interval> &, const std::vector<Interval> &);

bool is_rectangle(const std::vector<std::vector<long double>> &mat) {
  if (mat.empty())
    return true;

  for (size_t i = 0; i < mat.size(); i++)
    if (mat[0].size() != mat[i].size())
      return false;

  return true;
}

// if mat is not empty or not rectangle, returns pair of {num_rows, num_cols};
std::pair<int, int> get_mat_dimensions(const std::vector<std::vector<long double>> &mat) {
  if (!is_rectangle(mat) || mat.empty()) {
    logger.error("Matrix is not rectangle. Can't calculate dimensions.");
    exit(1);
  }

  return {mat.size(), mat[0].size()};
}

std::vector<std::vector<long double>> matrix_multiply(const std::vector<std::vector<long double>> &mat1,
                                                      const std::vector<std::vector<long double>> &mat2) {
  auto mat1_dim = get_mat_dimensions(mat1), mat2_dim = get_mat_dimensions(mat2);
  if (mat1_dim.second != mat2_dim.first) {
    logger.error("Matrix dimensions do not match for multiplication.");
    exit(1);
  }

  std::vector<std::vector<long double>> result(mat1_dim.first, std::vector<long double>(mat2_dim.second));

  for (int i = 0; i < mat1_dim.first; ++i) {
    for (int j = 0; j < mat2_dim.second; ++j) {
      for (int k = 0; k < mat1_dim.second; ++k) {
        result[i][j] += mat1[i][k] * mat2[k][j];
      }
    }
  }

  return result;
}

std::array<std::array<long double, 2>, 2> matrix_multiply(const std::array<std::array<long double, 2>, 2> &mat1,
                                                          const std::array<std::array<long double, 2>, 2> &mat2) {
  std::array<std::array<long double, 2>, 2> result{};

  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      for (int k : {0, 1}) {
        result[i][j] += mat1[i][k] * mat2[k][j];
      }
    }
  }

  return result;
}

// checks if given matrix is a square
bool is_square(const std::vector<std::vector<long double>> &mat) {
  size_t rows = mat.size();
  for (size_t row = 0; row < rows; row++)
    if (mat[row].size() != rows)
      return false;
  return true;
}

std::vector<std::vector<long double>> binary_exponentiation(const std::vector<std::vector<long double>> &mat,
                                                            long long power) {
  if (!is_square(mat)) {
    logger.error("Matrix must be square for exponentiation.");
    exit(1);
  }

  std::vector<std::vector<long double>> result(mat.size());
  for (size_t i = 0; i < result.size(); i++) {
    result[i].resize(mat.size());
    result[i][i] = 1;
  }

  std::vector<std::vector<long double>> base = mat;

  while (power > 0) {
    if (power % 2 == 1)
      result = matrix_multiply(result, base);
    base = matrix_multiply(base, base);
    power /= 2;
  }

  return result;
}

std::array<std::array<long double, 2>, 2> binary_exponentiation(const std::array<std::array<long double, 2>, 2> &mat,
                                                                long long power) {
  std::array<std::array<long double, 2>, 2> result{};
  for (size_t i = 0; i < result.size(); i++) {
    result[i][i] = 1;
  }

  std::array<std::array<long double, 2>, 2> base = mat;

  while (power > 0) {
    if (power % 2 == 1)
      result = matrix_multiply(result, base);
    base = matrix_multiply(base, base);
    power /= 2;
  }

  return result;
}

long double logsumexp(const std::vector<long double> &values) {
  if (values.empty()) {
    return 0.0L;
  }

  long double ld_inf = std::numeric_limits<long double>::infinity();

  long double max_value = values[0];
  for (size_t i = 1; i < values.size(); i++) {
    if (values[i] - max_value > 1e-50) {
      max_value = values[i];
    }
  }

  if (max_value == -ld_inf) {
    return max_value;
  }

  long double sum = 0.0L;
  for (long double value : values)
    sum += exp(value - max_value);

  return max_value + log(sum);
}

std::vector<long double> joint_logprobs(const std::vector<std::vector<long double>> &probs_by_chr) {
  if (probs_by_chr.size() == 0) {
    logger.error("p-values should have at least one level!.");
    exit(1);
  }

  for (size_t idx = 0; idx < probs_by_chr.size(); idx++) {
    if (probs_by_chr[idx].empty()) {
      logger.error("Layers should be non-empty!");
      exit(1);
    }
  }

  if (probs_by_chr.size() == 1) {
    return probs_by_chr[0];
  }

  long long max_k = 0;
  for (size_t idx = 0; idx < probs_by_chr.size(); idx++)
    max_k += probs_by_chr[idx].size() - 1;

  const long double ld_inf = std::numeric_limits<long double>::infinity();

  std::vector<long double> prev_row(max_k + 1, -ld_inf);
  for (size_t pos = 0; pos < probs_by_chr[0].size(); pos++)
    prev_row[pos] = probs_by_chr[0][pos];

  std::vector<long double> next_row(max_k + 1, -ld_inf);

  std::vector<long double> accum;
  for (std::vector<long double> level :
       std::vector<std::vector<long double>>(probs_by_chr.begin() + 1, probs_by_chr.end())) {
    for (long long k = 0; k <= max_k; k++) {
      long long accum_size = std::min(k + 1, (long long)level.size());
      accum.resize(accum_size);
      for (int j = 0; j < accum_size; j++)
        accum[j] = level[j] + prev_row[k - j];
      next_row[k] = logsumexp(accum);
      accum.clear();
    }
    prev_row = next_row;
  }

  return next_row;
}

MultiProbs joint_logprobs(const MultiProbs &probs1, const MultiProbs &probs2) {
  if (probs1.empty() || probs2.empty()) {
    logger.error("multi probs should not be empty");
    exit(1);
  }

  if (probs1.size() != 2 || probs2.size() != 2) {
    logger.error("invalid multi probs dimensions, 2x2 necessary");
    exit(1);
  }

  for (int i : {0, 1}) {
    if (probs1[i].size() != 2 || probs2[i].size() != 2) {
      logger.error("invalid multi probs dimensions, 2x2 necessary");
      exit(1);
    }
    for (int j : {0, 1}) {
      if (probs1[i][j].empty() || probs2[i][j].empty()) {
        logger.error("probs in multi probs should not be empty");
        exit(1);
      }
    }
  }

  MultiProbs res{};

  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      std::vector<long double> midpoint_0 = joint_logprobs({probs1[i][0], probs2[0][j]}),
                               midpoint_1 = joint_logprobs({probs1[i][1], probs2[1][j]});

      if (midpoint_0.size() != midpoint_1.size()) {
        logger.error("combined probs with different intermediary states do not "
                     "have the same lengths. this should not happen :D");
        exit(1);
      }

      std::vector<long double> combined(midpoint_0.size());
      for (size_t idx = 0; idx < midpoint_0.size(); idx++) {
        combined[idx] = logsumexp({midpoint_0[idx], midpoint_1[idx]});
      }
      res[i][j] = combined;
    }
  }

  return res;
}

// calculate joint p-value for a given `overlap_count`.
// `p_values_by_level` should contain log-values
// if significance is enrichment calculates right tail, otherwise calculates left tail
// combined is not allowed
long double calculate_joint_pvalue(const std::vector<std::vector<long double>> &probs_by_chr, long long overlap_count,
                                   Significance significance) {
  if (significance == Significance::COMBINED) {
    logger.error("Significance cannot be combined when calculating joint p-value.");
    exit(1);
  }

  if (overlap_count < 0 || probs_by_chr.empty() || (probs_by_chr.size() == 1 && probs_by_chr[0].empty()))
    return 1;

  std::vector<long double> logprobs = joint_logprobs(probs_by_chr);
  bool is_enrichment = significance == Significance::ENRICHMENT;
  if (is_enrichment && overlap_count >= (long long)logprobs.size())
    return 0;

  std::vector<long double> trimmed_probs;
  size_t start_idx = is_enrichment ? overlap_count : 0, end_idx = is_enrichment ? logprobs.size() : overlap_count + 1;
  for (size_t idx = start_idx; idx < end_idx; idx++) {
    trimmed_probs.push_back(logprobs[idx]);
  }
  long double result = exp(logsumexp(trimmed_probs));
  return result;
}

std::vector<std::vector<long double>> vector_to_2d_matrix(const std::vector<long double> &vec) {
  return std::vector<std::vector<long double>>{vec};
}

bool same_dimensions(const std::vector<std::vector<long double>> &mat1,

                     const std::vector<std::vector<long double>> &mat2) {

  auto mat1_dim = get_mat_dimensions(mat1), mat2_dim = get_mat_dimensions(mat2);
  return mat1_dim == mat2_dim;
}

std::vector<std::vector<long double>> add_matrices(const std::vector<std::vector<long double>> &mat1,
                                                   const std::vector<std::vector<long double>> &mat2) {
  if (!same_dimensions(mat1, mat2)) {
    logger.error("Can't add matrices. They don't have the same dimensions.");
    exit(1);
  }

  auto dim = get_mat_dimensions(mat1);
  std::vector<std::vector<long double>> result(dim.first, std::vector<long double>(dim.second));
  for (int i = 0; i < dim.first; i++)
    for (int j = 0; j < dim.second; j++)
      result[i][j] = mat1[i][j] + mat2[i][j];

  return result;
}

std::array<std::array<long double, 2>, 2> add_matrices(const std::array<std::array<long double, 2>, 2> &mat1,
                                                       const std::array<std::array<long double, 2>, 2> &mat2) {
  std::array<std::array<long double, 2>, 2> result{};
  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      result[i][j] = mat1[i][j] + mat2[i][j];
    }
  }

  return result;
}

std::vector<std::vector<long double>> subtract_matrices(const std::vector<std::vector<long double>> &mat1,
                                                        const std::vector<std::vector<long double>> &mat2) {

  if (!same_dimensions(mat1, mat2)) {
    logger.error("Can't subtract matrices. They don't have the same dimensions.");
    exit(1);
  }

  auto dim = get_mat_dimensions(mat1);
  std::vector<std::vector<long double>> result(dim.first, std::vector<long double>(dim.second));
  for (int i = 0; i < dim.first; i++)
    for (int j = 0; j < dim.second; j++)
      result[i][j] = mat1[i][j] - mat2[i][j];

  return result;
}

std::array<std::array<long double, 2>, 2> subtract_matrices(const std::array<std::array<long double, 2>, 2> &mat1,
                                                            const std::array<std::array<long double, 2>, 2> &mat2) {
  std::array<std::array<long double, 2>, 2> result{};
  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      result[i][j] = mat1[i][j] - mat2[i][j];
    }
  }

  return result;
}

std::vector<long double> matrix_to_vector(const std::vector<std::vector<long double>> &mat) {
  if (mat.size() != 1) {
    logger.error("Invalid matrix dimensions for converting into vector. Only "
                 "possible for matrices with dimensions 1xN");
    exit(1);
  }

  return mat[0];
}

std::string to_string(const std::vector<std::vector<long double>> &matrix) {
  std::ostringstream oss;
  oss << "[";

  for (size_t i = 0; i < matrix.size(); ++i) {
    if (i > 0)
      oss << ", ";
    oss << "[";

    for (size_t j = 0; j < matrix[i].size(); ++j) {
      if (j > 0)
        oss << ", ";
      oss << matrix[i][j];
    }

    oss << "]";
  }

  oss << "]";
  return oss.str();
}

std::string to_string(const std::vector<long double> &vec) {
  std::ostringstream oss;
  oss << "[";
  oss << std::fixed << std::setprecision(12);

  for (size_t i = 0; i < vec.size(); ++i) {
    if (i > 0)
      oss << ", ";
    oss << vec[i];
  }

  oss << "]";
  return oss.str();
}

std::string to_string(const std::vector<Interval> &vec) {
  std::ostringstream oss;
  oss << "[";

  for (size_t i = 0; i < vec.size(); ++i) {
    if (i > 0)
      oss << ", ";
    oss << vec[i];
  }

  oss << "]";
  return oss.str();
}

std::string to_string(const MultiProbs &multi_probs) {
  std::ostringstream oss;

  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      oss << "(" << i << "," << j << "): " << to_string(multi_probs[i][j]) << "\n";
    }
  }

  return oss.str();
}

std::string to_string(const long double &val) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(30) << val;
  return oss.str();
}

Interval slice_interval_by_window(const Interval &window, const Interval &interval) {
  return {interval.chr_name, std::max(window.begin, interval.begin), std::min(window.end, interval.end)};
}

bool are_intervals_non_overlapping(const std::vector<Interval> &intervals) {
  std::vector<Interval> sortedIntervals = intervals;
  std::sort(sortedIntervals.begin(), sortedIntervals.end());
  for (size_t idx = 1; idx < sortedIntervals.size(); idx++) {
    if (sortedIntervals[idx].begin < sortedIntervals[idx - 1].end) {
      return false;
    }
  }

  return true;
}

WindowSectionSplitResult split_windows_into_non_overlapping_sections(const std::vector<Interval> &windows,
                                                                     const std::vector<Interval> &ref_intervals,
                                                                     const std::vector<Interval> &query_intervals) {
  if (windows.empty()) {
    return {};
  }

  std::string chr_name = windows[0].chr_name;

  const int REF_INTERVAL = 0;
  const int QUERY_INTERVAL = 1;
  const int WINDOW = 2;
  struct Event {
    long long pos;
    bool end;
    size_t idx;
    int type;
  };
  std::vector<Event> events(2 * windows.size() + 2 * ref_intervals.size() + 2 * query_intervals.size());
  for (size_t idx = 0; idx < windows.size(); idx++) {
    Interval window = windows[idx];
    events[idx << 1] = {window.begin, 0, idx, WINDOW};
    events[(idx << 1) | 1] = {window.end, 1, idx, WINDOW};
  }

  size_t buf = 2 * windows.size();
  for (size_t idx = 0; idx < ref_intervals.size(); idx++) {
    Interval ref_interval = ref_intervals[idx];
    events[buf + (idx << 1)] = {ref_interval.begin, 0, idx, REF_INTERVAL};
    events[buf + ((idx << 1) | 1)] = {ref_interval.end, 1, idx, REF_INTERVAL};
  }

  buf += 2 * ref_intervals.size();
  for (size_t idx = 0; idx < query_intervals.size(); idx++) {
    Interval query_interval = query_intervals[idx];
    events[buf + (idx << 1)] = {query_interval.begin, 0, idx, QUERY_INTERVAL};
    events[buf + ((idx << 1) | 1)] = {query_interval.end, 1, idx, QUERY_INTERVAL};
  }

  std::sort(events.begin(), events.end(), [](const Event &e1, const Event &e2) {
    if (e1.pos == e2.pos) {
      if (e1.type == e2.type)
        return e1.type != WINDOW ? e1.end < e2.end : e1.end > e2.end;
      if (e1.type == WINDOW || e2.type == WINDOW) {
        return e1.end == e2.end ? e1.type < e2.type : e1.end > e2.end;
      }
    }
    return e1.pos < e2.pos;
  });

  std::vector<Interval> spans(windows.size());
  long long last_pos = -1, section_start = -1;

  std::vector<Section> sections;
  Interval currently_opened_ref_interval, currently_opened_query_interval;
  Interval start_ref_int, start_query_int;
  long long opened = 0;
  for (Event event : events) {
    if (event.type == REF_INTERVAL) {
      Interval current_ref_interval = ref_intervals[event.idx];
      if (!event.end) {
        currently_opened_ref_interval = current_ref_interval;
      } else {
        currently_opened_ref_interval = Interval("", -1, -1);
      }
    } else if (event.type == QUERY_INTERVAL) {
      Interval current_query_interval = query_intervals[event.idx];
      if (!event.end) {
        currently_opened_query_interval = current_query_interval;
      } else {
        currently_opened_query_interval = Interval("", -1, -1);
      }
    } else {
      if (opened > 0 && section_start != -1 && event.pos != last_pos) {
        long long section_end = event.pos;
        long long ref_interval_begin = start_ref_int.length() != 0 ? start_ref_int.get_begin()
                                       : currently_opened_ref_interval.length() != 0
                                           ? currently_opened_ref_interval.get_begin()
                                           : section_start;
        long long ref_interval_end =
            currently_opened_ref_interval.length() != 0 ? currently_opened_ref_interval.get_end() : section_end;
        long long query_interval_begin = start_query_int.length() != 0 ? start_query_int.get_begin()
                                         : currently_opened_query_interval.length() != 0
                                             ? currently_opened_query_interval.get_begin()
                                             : section_start;
        long long query_interval_end =
            currently_opened_query_interval.length() != 0 ? currently_opened_query_interval.get_end() : section_end;
        sections.push_back(Section(chr_name, section_start, section_end, ref_interval_begin < section_start,
                                   section_end < ref_interval_end, query_interval_begin < section_start,
                                   section_end < query_interval_end));
      }
      if (event.pos != last_pos) {
        section_start = event.pos;
      }
      last_pos = event.pos;

      if (!event.end) {
        opened++;
        spans[event.idx].begin = sections.size();
      } else {
        opened--;
        spans[event.idx].chr_name = chr_name;
        spans[event.idx].end = sections.size();
      }

      start_ref_int = currently_opened_ref_interval;
      start_query_int = currently_opened_query_interval;
    }
  }

  return WindowSectionSplitResult(sections, spans);
}

long double log_multiply(long double log_x, long double log_y) {
  const long double ld_inf = std::numeric_limits<long double>::max();
  if (log_x == -ld_inf)
    return log_x;
  if (log_y == -ld_inf)
    return log_y;
  return log_x + log_y;
}

std::vector<long double> merge_multi_probs(MultiProbs probs, const MarkovChain &markov_chain) {
  if (probs.size() != 2 || probs[0].size() != 2 || probs[1].size() != 2) {
    logger.error("invalid multiprobs or stationary_distribution "
                 "dimensions/size for merging into single");
    exit(1);
  }

  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      if (probs[0][0].size() != probs[i][j].size()) {
        logger.error("multi probs probs are not the same length");
        exit(1);
      }
    }
  }

  size_t k = probs[0][0].size();
  std::vector<long double> res(k);

  for (size_t i : {0, 1}) {
    for (size_t idx = 0; idx < k; idx++) {
      probs[i][0][idx] = logsumexp({probs[i][0][idx], probs[i][1][idx]});
    }
  }

  StationaryDistribution stationary_distribution = markov_chain.get_stationary_distribution();

  for (size_t idx = 0; idx < k; idx++) {
    res[idx] = logsumexp({log_multiply(logl(stationary_distribution[0]), probs[0][0][idx]),
                          log_multiply(logl(stationary_distribution[1]), probs[1][0][idx])});
  }

  return res;
}

void print_multiprobs(const MultiProbs &probs) {
  if (probs.size() != 2 || probs[0].size() != 2 || probs[1].size() != 2) {
    std::cout << "invalid dimensions of multiprobs for printing\n";
    return;
  }

  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      std::cout << "(" << i << "," << j << "): " << probs[i][j] << "\n";
    }
  }
}

Section join_sections(const Section &section1, const Section &section2, const MarkovChain &markov_chain) {
  SectionProbs probs1 = section1.get_probs(), probs2 = section2.get_probs();
  std::vector<Interval> ref_ints1 = section1.get_ref_intervals(), ref_ints2 = section2.get_ref_intervals();
  std::vector<Interval> query_ints1 = section1.get_query_intervals(), query_ints2 = section2.get_query_intervals();

  bool ref_overflows = section1.get_last_ref_interval_intersected() && section2.get_first_ref_interval_intersected();
  bool query_overflows =
      section1.get_last_query_interval_intersected() && section2.get_first_query_interval_intersected();
  MultiProbs middle_probs;
  std::vector<Interval> ref_intervals, query_intervals;
  if (ref_overflows) {
    Interval last_ref_section1 = ref_ints1.back();
    Interval first_ref_section2 = ref_ints2.front();
    ref_intervals = {
        Interval(last_ref_section1.get_chr_name(), last_ref_section1.get_begin(), first_ref_section2.get_end())};
    long long new_start = last_ref_section1.get_begin();
    middle_probs =
        Model::eval_probs_single_chr_direct_new(ref_intervals, new_start, first_ref_section2.get_end(), markov_chain);
  }

  if (query_overflows) {
    Interval last_query_section1 = query_ints1.back();
    Interval first_query_section2 = query_ints2.front();
    query_intervals = {
        Interval(last_query_section1.get_chr_name(), last_query_section1.get_begin(), first_query_section2.get_end())};
  }

  bool should_calculate_middle_probs = !ref_ints1.empty() && !ref_ints2.empty() &&
                                       ref_ints1.back().get_begin() != section1.get_begin() &&
                                       ref_ints2.front().get_end() != section2.get_end();
  bool last_fills_first = !ref_ints1.empty() && ref_ints1.back().get_begin() == section1.get_begin();
  bool first_fills_second = !ref_ints2.empty() && ref_ints2.front().get_end() == section2.get_end();

  MultiProbs new_normal = joint_logprobs(probs1.get_normal(), probs2.get_normal());
  if (ref_overflows) {
    new_normal = joint_logprobs(joint_logprobs(probs1.get_except_last(), middle_probs), probs2.get_except_first());
  }

  MultiProbs new_except_first = joint_logprobs(probs1.get_except_first(), probs2.get_normal());
  if (ref_overflows) {
    new_except_first = should_calculate_middle_probs
                           ? joint_logprobs(joint_logprobs(probs1.get_except_first_and_last(), middle_probs),
                                            probs2.get_except_first())
                       : last_fills_first ? probs2.get_except_first()
                                          : joint_logprobs(probs1.get_except_first_and_last(), middle_probs);
  }

  MultiProbs new_except_last = joint_logprobs(probs1.get_normal(), probs2.get_except_last());
  if (ref_overflows) {
    new_except_last =
        should_calculate_middle_probs
            ? joint_logprobs(joint_logprobs(probs1.get_except_last(), middle_probs), probs2.get_except_first_and_last())
        : first_fills_second ? probs1.get_except_last()
                             : joint_logprobs(middle_probs, probs2.get_except_first_and_last());
  }

  MultiProbs new_except_first_and_last = joint_logprobs(probs1.get_except_first(), probs2.get_except_last());
  if (ref_overflows) {
    new_except_first_and_last = should_calculate_middle_probs
                                    ? joint_logprobs(joint_logprobs(probs1.get_except_first_and_last(), middle_probs),
                                                     probs2.get_except_first_and_last())
                                : last_fills_first ? probs2.get_except_first_and_last()
                                                   : probs1.get_except_first_and_last();
  }

  long long new_overlap_count = section1.get_overlap_count() + section2.get_overlap_count();
  if (ref_overflows && !query_ints1.empty() && !ref_ints1.empty() &&
      query_ints1.back().get_end() > ref_ints1.back().get_begin() && !query_ints2.empty() && !ref_ints2.empty() &&
      query_ints2.front().get_begin() < ref_ints2.front().get_end()) {
    new_overlap_count--;
  }

  std::vector<Interval> new_ref_intervals = ref_ints1;
  std::vector<Interval> new_query_intervals = query_ints1;
  if (!ref_overflows) {
    new_ref_intervals.insert(new_ref_intervals.end(), ref_ints2.begin(), ref_ints2.end());
  } else {
    new_ref_intervals.pop_back();
    new_ref_intervals.push_back(ref_intervals[0]);
    new_ref_intervals.insert(new_ref_intervals.end(), ref_ints2.begin() + (!ref_ints2.empty()), ref_ints2.end());
  }

  if (!query_overflows) {
    new_query_intervals.insert(new_query_intervals.end(), query_ints2.begin(), query_ints2.end());
  } else {
    new_query_intervals.pop_back();
    new_query_intervals.push_back(query_intervals[0]);
    new_query_intervals.insert(new_query_intervals.end(), query_ints2.begin() + (!query_ints2.empty()),
                               query_ints2.end());
  }

  Section merged_section(section1.get_chr_name(), section1.get_begin(), section2.get_end(),
                         section1.get_first_ref_interval_intersected(), section2.get_last_ref_interval_intersected(),
                         section1.get_first_query_interval_intersected(),
                         section2.get_last_query_interval_intersected(), new_ref_intervals, new_query_intervals);
  SectionProbs merged_probs(new_normal, new_except_first, new_except_last, new_except_first_and_last);
  merged_section.set_probs(merged_probs);
  merged_section.set_overlap_count(new_overlap_count);

  return merged_section;
}

bool compare_logprobs_vectors(const std::vector<long double> &a, const std::vector<long double> &b,
                              long double epsilon) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(), [epsilon](long double x, long double y) {
           return std::abs(exp(x) - exp(y)) < epsilon;
         });
}

bool compare_multiprobs(const MultiProbs &a, const MultiProbs &b, long double epsilon) {
  for (int i : {0, 1}) {
    for (int j : {0, 1}) {
      if (!compare_logprobs_vectors(a[i][j], b[i][j]))
        return false;
    }
  }
  return true;
}

std::vector<Interval> split_intervals_into_ones(const std::vector<Interval> &intervals) {
  std::vector<Interval> new_intervals;

  for (Interval interval : intervals) {
    for (long long pos = interval.get_begin(); pos < interval.get_end(); pos++) {
      new_intervals.push_back(Interval(interval.get_chr_name(), pos, pos + 1));
    }
  }

  return new_intervals;
}

Section join_sections_new(const Section &section1, const Section &section2, const MarkovChain &markov_chain) {
  SectionProbs probs1 = section1.get_probs(), probs2 = section2.get_probs();
  std::vector<Interval> ref_ints1 = section1.get_ref_intervals(), ref_ints2 = section2.get_ref_intervals();
  std::vector<Interval> query_ints1 = section1.get_query_intervals(), query_ints2 = section2.get_query_intervals();

  bool ref_overflows = !ref_ints1.empty() && section1.get_last_ref_interval_intersected() && !ref_ints2.empty() &&
                       section2.get_first_ref_interval_intersected();
  bool query_overflows = !query_ints1.empty() && section1.get_last_query_interval_intersected() &&
                         !query_ints2.empty() && section2.get_first_query_interval_intersected();
  MultiProbs middle_probs;
  std::vector<Interval> ref_intervals, query_intervals;
  if (ref_overflows) {
    Interval last_ref_section1 = ref_ints1.back();
    Interval first_ref_section2 = ref_ints2.front();
    ref_intervals = {
        Interval(last_ref_section1.get_chr_name(), last_ref_section1.get_begin(), first_ref_section2.get_end())};
    middle_probs = Model::eval_probs_single_chr_direct_new(ref_intervals, last_ref_section1.get_begin(),
                                                           first_ref_section2.get_end(), markov_chain);
  }

  if (query_overflows) {
    Interval last_query_section1 = query_ints1.back();
    Interval first_query_section2 = query_ints2.front();
    query_intervals = {
        Interval(last_query_section1.get_chr_name(), last_query_section1.get_begin(), first_query_section2.get_end())};
  }

  bool should_calculate_middle_probs = !ref_ints1.empty() && !ref_ints2.empty() &&
                                       ref_ints1.back().get_begin() != section1.get_begin() &&
                                       ref_ints2.front().get_end() != section2.get_end();
  bool last_fills_first = !ref_ints1.empty() && ref_ints1.back().get_begin() == section1.get_begin();

  MultiProbs new_probs = joint_logprobs(probs1.get_except_first_and_last(), probs2.get_except_first_and_last());
  if (ref_overflows) {
    if (should_calculate_middle_probs) {
      new_probs = joint_logprobs(joint_logprobs(probs1.get_except_first_and_last(), middle_probs),
                                 probs2.get_except_first_and_last());
    } else {
      if (last_fills_first) {
        new_probs = probs2.get_except_first_and_last();
        if (!section1.get_first_ref_interval_intersected())
          new_probs = joint_logprobs(middle_probs, new_probs);
      } else {
        new_probs = probs1.get_except_first_and_last();
        if (!section2.get_last_ref_interval_intersected())
          new_probs = joint_logprobs(new_probs, middle_probs);
      }
    }
  }

  long long new_overlap_count = section1.get_overlap_count() + section2.get_overlap_count();
  if (ref_overflows && !query_ints1.empty() && !ref_ints1.empty() &&
      query_ints1.back().get_end() > ref_ints1.back().get_begin() && !query_ints2.empty() && !ref_ints2.empty() &&
      query_ints2.front().get_begin() < ref_ints2.front().get_end()) {
    new_overlap_count--;
  }

  std::vector<Interval> new_ref_intervals = ref_ints1;
  std::vector<Interval> new_query_intervals = query_ints1;
  if (!ref_overflows) {
    new_ref_intervals.insert(new_ref_intervals.end(), ref_ints2.begin(), ref_ints2.end());
  } else {
    new_ref_intervals.pop_back();
    new_ref_intervals.push_back(ref_intervals[0]);
    new_ref_intervals.insert(new_ref_intervals.end(), ref_ints2.begin() + (!ref_ints2.empty()), ref_ints2.end());
  }

  if (!query_overflows) {
    new_query_intervals.insert(new_query_intervals.end(), query_ints2.begin(), query_ints2.end());
  } else {
    new_query_intervals.pop_back();
    new_query_intervals.push_back(query_intervals[0]);
    new_query_intervals.insert(new_query_intervals.end(), query_ints2.begin() + (!query_ints2.empty()),
                               query_ints2.end());
  }

  Section merged_section(section1.get_chr_name(), section1.get_begin(), section2.get_end(),
                         section1.get_first_ref_interval_intersected(), section2.get_last_ref_interval_intersected(),
                         section1.get_first_query_interval_intersected(),
                         section2.get_last_query_interval_intersected(), new_ref_intervals, new_query_intervals);
  SectionProbs merged_probs({}, {}, {}, new_probs);
  merged_section.set_probs(merged_probs);
  merged_section.set_overlap_count(new_overlap_count);

  return merged_section;
}

template <class T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec) {
  out << "[";
  for (size_t i = 0; i < vec.size(); i++) {
    if (i)
      out << " ";
    out << vec[i];
  }
  out << "]";
  return out;
}

Section join_sections_segtree(const Section &section1, const Section &section2, const MarkovChain &markov_chain) {
  if (section1.length() == 0)
    return section2;
  if (section2.length() == 0)
    return section1;
  return join_sections(section1, section2, markov_chain);
}

Section join_sections_new_segtree(const Section &section1, const Section &section2, const MarkovChain &markov_chain) {
  if (section1.length() == 0)
    return section2;
  if (section2.length() == 0)
    return section1;
  return join_sections_new(section1, section2, markov_chain);
}
