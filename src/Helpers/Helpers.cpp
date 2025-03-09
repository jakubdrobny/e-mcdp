#include "Helpers.hpp"
#include "../Interval/Interval.hpp"
#include "../Logger/Logger.hpp"
#include <algorithm>
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
                 str + " split by " + delimeter + ".");
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
                 ". Should be {chr_name} {begin} {end}.");
    exit(1);
  }

  return Interval(vals[0], std::stoll(vals[1]), std::stoll(vals[2]));
}

std::vector<Interval> load_intervals(const std::string &file_path,
                                     bool is_closed) {
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
std::vector<Interval>
load_windows(Args &args,
             std::unordered_map<std::string, long long> &chr_sizes) {
  std::vector<Interval> windows;

  if (args.windows_source == "basic" || args.windows_source == "dense") {
    std::vector<std::pair<std::string, long long>> chr_sizes_vec(
        chr_sizes.begin(), chr_sizes.end());
    std::sort(chr_sizes_vec.begin(), chr_sizes_vec.end());
    for (std::pair<std::string, long long> chr : chr_sizes) {
      std::string chr_name = chr.first;
      long long chr_size = chr.second;
      long long l = 0, r = args.windows_size;
      while (l < chr_size) {
        windows.push_back({chr_name, l, r});
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
      logger.error("Size for chromosome " + chr_name +
                   " specified more than once.");
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
    for (; idx < (int)intervals.size() && intervals[idx].chr_name <= chr_name;
         idx++)
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

std::vector<std::vector<long double>>
get_base_transition_matrix(long long chr_size,
                           std::vector<Interval> &query_intervals) {
  std::vector<std::vector<long double>> t(2, std::vector<long double>(2));

  long double L = chr_size;
  long double len_Q = query_intervals.size();
  long double weight_Q = 0;
  for (Interval interval : query_intervals)
    weight_Q += interval.end - interval.begin;

  t[0][1] = (len_Q) / (L - weight_Q - 1);
  t[0][0] = 1 - t[0][1];

  t[1][0] = (len_Q) / (weight_Q);
  t[1][1] = 1 - t[1][0];

  return t;
}

// returns T and T_mod (its just T, with zeros in second col)
std::pair<std::vector<std::vector<long double>>,
          std::vector<std::vector<long double>>>
get_transition_matrices(long long chr_size,
                        std::vector<Interval> &query_intervals) {
  if (query_intervals.empty()) {
    logger.error("Query intervals should not be empty.");
    exit(1);
  }

  std::vector<std::vector<long double>> t =
      get_base_transition_matrix(chr_size, query_intervals);
  std::vector<std::vector<long double>> d(2, std::vector<long double>(2));
  d[0][0] = t[0][0];
  d[1][0] = t[1][0];

  return {t, d};
}

template void extend<Interval>(std::vector<Interval> &,
                               const std::vector<Interval> &);

bool is_rectangle(const std::vector<std::vector<long double>> &mat) {
  if (mat.empty())
    return true;

  for (size_t i = 0; i < mat.size(); i++)
    if (mat[0].size() != mat[i].size())
      return false;

  return true;
}

// if mat is not empty or not rectangle, returns pair of {num_rows, num_cols};
std::pair<int, int>
get_mat_dimensions(const std::vector<std::vector<long double>> &mat) {
  if (!is_rectangle(mat) || mat.empty()) {
    logger.error("Matrix is not rectangle. Can't calculate dimensions.");
    exit(1);
  }

  return {mat.size(), mat[0].size()};
}

std::vector<std::vector<long double>>
matrix_multiply(const std::vector<std::vector<long double>> &mat1,
                const std::vector<std::vector<long double>> &mat2) {
  auto mat1_dim = get_mat_dimensions(mat1), mat2_dim = get_mat_dimensions(mat2);
  if (mat1_dim.second != mat2_dim.first) {
    logger.error("Matrix dimensions do not match for multiplication.");
    exit(1);
  }

  std::vector<std::vector<long double>> result(
      mat1_dim.first, std::vector<long double>(mat2_dim.second));

  for (int i = 0; i < mat1_dim.first; ++i) {
    for (int j = 0; j < mat2_dim.second; ++j) {
      for (int k = 0; k < mat1_dim.second; ++k) {
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

std::vector<std::vector<long double>>
binary_exponentiation(const std::vector<std::vector<long double>> &mat,
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

long double logsumexp(const std::vector<long double> &values) {
  const long double ld_inf = std::numeric_limits<long double>::infinity();
  if (values.size() == 0)
    return -ld_inf;

  long double max_value = 0.0L;
  for (size_t i = 0; i < values.size(); i++)
    max_value = std::max(max_value, values[i]);

  long double sum = 0.0L;
  for (long double value : values)
    sum += exp(value - max_value);

  return max_value + log(sum);
}

std::vector<long double>
joint_logprobs(const std::vector<std::vector<long double>> &probs_by_chr) {
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
  for (std::vector<long double> level : std::vector<std::vector<long double>>(
           probs_by_chr.begin() + 1, probs_by_chr.end())) {
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

// calculate joint p-value for a given `overlap_count`.
// `p_values_by_level` should contain log-values
long double calculate_joint_pvalue(
    const std::vector<std::vector<long double>> &probs_by_chr,
    long long overlap_count) {
  if (overlap_count < 0)
    return 1;

  std::vector<long double> logprobs = joint_logprobs(probs_by_chr);
  if (overlap_count >= (long long)logprobs.size())
    return 0;

  long double result = exp(logsumexp(std::vector<long double>(
      logprobs.begin() + overlap_count, logprobs.end())));
  return result;
}

std::vector<std::vector<long double>>
vector_to_2d_matrix(const std::vector<long double> &vec) {
  return std::vector<std::vector<long double>>{vec};
}

bool same_dimensions(const std::vector<std::vector<long double>> &mat1,

                     const std::vector<std::vector<long double>> &mat2) {

  auto mat1_dim = get_mat_dimensions(mat1), mat2_dim = get_mat_dimensions(mat2);
  return mat1_dim == mat2_dim;
}

std::vector<std::vector<long double>>
add_matrices(const std::vector<std::vector<long double>> &mat1,
             const std::vector<std::vector<long double>> &mat2) {
  if (!same_dimensions(mat1, mat2)) {
    logger.error("Can't add matrices. They don't have the same dimensions.");
    exit(1);
  }

  auto dim = get_mat_dimensions(mat1);
  std::vector<std::vector<long double>> result(
      dim.first, std::vector<long double>(dim.second));
  for (int i = 0; i < dim.first; i++)
    for (int j = 0; j < dim.second; j++)
      result[i][j] = mat1[i][j] + mat2[i][j];

  return result;
}

std::vector<std::vector<long double>>
subtract_matrices(const std::vector<std::vector<long double>> &mat1,
                  const std::vector<std::vector<long double>> &mat2) {

  if (!same_dimensions(mat1, mat2)) {
    logger.error(
        "Can't subtract matrices. They don't have the same dimensions.");
    exit(1);
  }

  auto dim = get_mat_dimensions(mat1);
  std::vector<std::vector<long double>> result(
      dim.first, std::vector<long double>(dim.second));
  for (int i = 0; i < dim.first; i++)
    for (int j = 0; j < dim.second; j++)
      result[i][j] = mat1[i][j] - mat2[i][j];

  return result;
}

std::vector<long double>
matrix_to_vector(const std::vector<std::vector<long double>> &mat) {
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

  for (size_t i = 0; i < vec.size(); ++i) {
    if (i > 0)
      oss << ", ";
    oss << vec[i];
  }

  oss << "]";
  return oss.str();
}

std::string to_string(const long double &val) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(30) << val;
  return oss.str();
}

// expects transition matrix 2x2
std::vector<long double>
get_stationary_distribution(const std::vector<std::vector<long double>> &mat) {
  if (mat.size() != 2 || mat[0].size() != 2 || mat[1].size() != 2) {
    logger.error("Invalid transition matrix dimensions. Can't calculate "
                 "stationary distribution.");
    exit(1);
  }

  // b = mat[0][1], a = 1 - b, c = mat[1][0], d = 1 - c
  // ax + cy = x and bx + dy = y and pi_0 + pi_1 = 1 should hold
  // solving for pi_0 and pi_1 we get:
  // we can derive that pi_0 = c/(b+c) and pi_1 = b/(b+c)
  long double b = mat[0][1], c = mat[1][0];
  long double denom = b + c;

  // is not irreducible
  if (std::abs(denom) < 1e-9) {
    logger.error(
        "Can't calculate stationary distribution. P[0][1] or P[1][0] are 0.");
    exit(1);
  }

  return {c / denom, b / denom};
}
