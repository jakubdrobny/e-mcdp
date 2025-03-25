#include "WindowResult.hpp"
#include "../Interval/Interval.hpp"
#include <vector>

WindowResult::WindowResult() {}

WindowResult::WindowResult(Interval window, long long overlap_count, std::vector<long double> probs)
    : window(window), overlap_count(overlap_count), probs(probs) {}

WindowResult::WindowResult(Interval window, long long overlap_count, MultiProbs multi_probs)
    : window(window), overlap_count(overlap_count), multi_probs(multi_probs) {}

Interval WindowResult::get_window() { return window; }
long long WindowResult::get_overlap_count() { return overlap_count; }
std::vector<long double> WindowResult::get_probs() { return probs; }
MultiProbs WindowResult::get_multi_probs() { return multi_probs; }
