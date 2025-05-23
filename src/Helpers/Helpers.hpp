#ifndef HELPERS_H
#define HELPERS_H

#include "../Args/Args.hpp"
#include "../Interval/Interval.hpp"
#include "../MarkovChain/MarkovChain.hpp"
#include "../Results/WindowResult.hpp"
#include "../Results/WindowSectionSplitResult.hpp"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using ChrSizesVector = std::vector<std::pair<std::string, long long>>;
using ChrSizesMap = std::unordered_map<std::string, long long>;

Interval parse_intervals_line(std::string line);

std::vector<Interval> load_intervals(const std::string &file_path, bool is_closed = false);

ChrSizesMap load_chr_sizes(const std::string &file_path);

std::unordered_set<std::string>
load_chr_names_from_chr_sizes(const std::unordered_map<std::string, long long> &chr_sizes);

std::vector<Interval> filter_intervals_by_chr_name(std::vector<Interval> intervals,
                                                   std::unordered_set<std::string> chr_names);

std::vector<Interval> merge_non_disjoint_intervals(std::vector<Interval> intervals);

std::vector<Interval> remove_empty_intervals(std::vector<Interval> intervals);

long long count_overlaps(std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals);

long long count_overlaps_single_chr(std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals);

std::vector<std::string> get_sorted_chr_names_from_intervals(std::vector<Interval> intervals);

ChrSizesVector chr_sizes_map_to_array(std::unordered_map<std::string, long long> &chr_sizes);

long double calculate_joint_pvalue(const std::vector<std::vector<long double>> &probs_by_chr, long long overlap_count,
                                   Significance significance);

std::vector<std::vector<long double>> get_base_transition_matrix(long long chr_size,
                                                                 const std::vector<Interval> &query_intervals);

std::pair<std::vector<std::vector<long double>>, std::vector<std::vector<long double>>>
get_transition_matrices(long long chr_size, const std::vector<Interval> &query_intervals);

template <typename T> void extend(std::vector<T> &self, const std::vector<T> &other);

std::vector<std::vector<long double>> matrix_multiply(const std::vector<std::vector<long double>> &mat1,
                                                      const std::vector<std::vector<long double>> &mat2);

std::array<std::array<long double, 2>, 2> matrix_multiply(const std::array<std::array<long double, 2>, 2> &mat1,
                                                          const std::array<std::array<long double, 2>, 2> &mat2);

std::vector<std::vector<long double>> binary_exponentiation(const std::vector<std::vector<long double>> &mat,
                                                            long long power);

std::array<std::array<long double, 2>, 2> binary_exponentiation(const std::array<std::array<long double, 2>, 2> &mat,
                                                                long long power);
std::vector<long double> joint_logprobs(const std::vector<std::vector<long double>> &probs_by_chr);

// joins two sets of log probs
MultiProbs joint_logprobs(const MultiProbs &probs1, const MultiProbs &probs2);

long double logsumexp(const std::vector<long double> &arr);

bool is_rectangle(const std::vector<std::vector<long double>> &mat);

std::pair<int, int> get_mat_dimensions(const std::vector<std::vector<long double>> &mat);

std::vector<std::vector<long double>> vector_to_2d_matrix(const std::vector<long double> &vec);

std::vector<std::vector<long double>> add_matrices(const std::vector<std::vector<long double>> &mat1,
                                                   const std::vector<std::vector<long double>> &mat2);

std::array<std::array<long double, 2>, 2> add_matrices(const std::array<std::array<long double, 2>, 2> &mat1,
                                                       const std::array<std::array<long double, 2>, 2> &mat2);

std::vector<std::vector<long double>> subtract_matrices(const std::vector<std::vector<long double>> &mat1,
                                                        const std::vector<std::vector<long double>> &mat2);

std::array<std::array<long double, 2>, 2> subtract_matrices(const std::array<std::array<long double, 2>, 2> &mat1,
                                                            const std::array<std::array<long double, 2>, 2> &mat2);

std::vector<long double> matrix_to_vector(const std::vector<std::vector<long double>> &mat);

std::string to_string(const std::vector<std::vector<long double>> &matrix);

std::string to_string(const std::vector<long double> &vec);

std::string to_string(const long double &val);

std::string to_string(const std::vector<Interval> &vec);

std::string to_string(const MultiProbs &multi_probs);

std::vector<Interval> load_windows(Args &args, std::unordered_map<std::string, long long> &chr_sizes);

std::vector<std::vector<Interval>> get_windows_intervals(const std::vector<Interval> &windows,
                                                         const std::vector<Interval> &intervals);

Interval slice_interval_by_window(const Interval &window, const Interval &interval);

bool are_intervals_non_overlapping(const std::vector<Interval> &intervals);

WindowSectionSplitResult split_windows_into_non_overlapping_sections(const std::vector<Interval> &windows,
                                                                     const std::vector<Interval> &ref_intervals,
                                                                     const std::vector<Interval> &query_intervals);

std::vector<long double> merge_multi_probs(MultiProbs probs, const MarkovChain &markov_chain);

void print_multiprobs(const MultiProbs &probs);

Section join_sections(const Section &section1, const Section &section2, const MarkovChain &markov_chain);
Section join_sections_new(const Section &section1, const Section &section2, const MarkovChain &markov_chain);

Section join_sections_segtree(const Section &section1, const Section &section2, const MarkovChain &markov_chain);
Section join_sections_new_segtree(const Section &section1, const Section &section2, const MarkovChain &markov_chain);

bool compare_logprobs_vectors(const std::vector<long double> &a, const std::vector<long double> &b,
                              long double epsilon = 1e-12);

bool compare_multiprobs(const MultiProbs &a, const MultiProbs &b, long double epsilon = 1e-12);

std::vector<Interval> split_intervals_into_ones(const std::vector<Interval> &intervals);

template <class T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &vec);

#endif
