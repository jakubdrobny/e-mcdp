#include "WindowResult.hpp"
#include "../Helpers/Helpers.hpp"
#include "../Interval/Interval.hpp"
#include <vector>

WindowResult::WindowResult() {}

WindowResult::WindowResult(Interval window, long long overlap_count, std::vector<long double> probs)
    : window(window), overlap_count(overlap_count), probs(probs) {}

WindowResult::WindowResult(Interval window, long long overlap_count, MultiProbs multi_probs)
    : window(window), overlap_count(overlap_count), multi_probs(multi_probs) {}

Interval WindowResult::get_window() const { return window; }

long long WindowResult::get_overlap_count() const { return overlap_count; }

std::vector<long double> WindowResult::get_probs() const { return probs; }

MultiProbs WindowResult::get_multi_probs() const { return multi_probs; }

bool WindowResult::operator==(const WindowResult &other) const {
  return this->window == other.window && compare_vectors_stl(this->probs, other.probs) &&
         this->overlap_count == other.overlap_count;
}

std::ostream &operator<<(std::ostream &os, const WindowResult &result) {
  os << "result:\n";
  os << "result.window: " << result.get_window() << "\n";
  os << "result.probs: " << to_string(result.get_probs()) << "\n";
  os << "result.overlap_count: " << result.get_overlap_count() << "\n";
  return os;
}
