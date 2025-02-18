#ifndef HELPERS_H
#define HELPERS_H

#include "../Interval/Interval.h"
#include "NumCpp.hpp"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using ChrSizesVector = std::vector<std::pair<std::string, long long>>;
using ChrSizesMap = std::unordered_map<std::string, long long>;

Interval parse_intervals_line(std::string line);

std::vector<Interval> load_intervals(const std::string &file_path,
                                     bool is_closed = false);

ChrSizesMap load_chr_sizes(const std::string &file_path);

std::unordered_set<std::string> load_chr_names_from_chr_sizes(
    const std::unordered_map<std::string, long long> &chr_sizes);

std::vector<Interval>
filter_intervals_by_chr_name(std::vector<Interval> intervals,
                             std::unordered_set<std::string> chr_names);

std::vector<Interval>
merge_non_disjoint_intervals(std::vector<Interval> intervals);

std::vector<Interval> remove_empty_intervals(std::vector<Interval> intervals);

long long count_overlaps(std::vector<Interval> ref_intervals,
                         std::vector<Interval> query_intervals);

std::vector<std::string>
get_sorted_chr_names_from_intervals(std::vector<Interval> intervals);

ChrSizesVector
chr_sizes_map_to_array(std::unordered_map<std::string, long long> &chr_sizes);

long double calculate_joint_pvalue(
    std::vector<std::vector<long double>> &probs_by_chromosome,
    long long overlap_count);

nc::NdArray<long double>
get_base_transition_matrix(long long chr_size,
                           std::vector<Interval> &query_intervals);

std::pair<nc::NdArray<long double>, nc::NdArray<long double>>
get_transition_matrices(long long chr_size,
                        std::vector<Interval> &query_intervals);

template <typename T>
void extend(std::vector<T> &self, const std::vector<T> &other);

nc::NdArray<long double> matrix_multiply(const nc::NdArray<long double> &mat1,
                                         const nc::NdArray<long double> &mat2);

nc::NdArray<long double>
binary_exponentiation(const nc::NdArray<long double> &mat, long long power);

nc::NdArray<long double>
matrix_row_to_2d_matrix(const nc::NdArray<long double> &mat, size_t row_index);
#endif
