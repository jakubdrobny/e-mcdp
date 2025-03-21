#ifndef WINDOWMODEL_H
#define WINDOWMODEL_H

#include "../Enums/Enums.hpp"
#include "../Helpers/Helpers.hpp"
#include "../Interval/Interval.hpp"
#include "../WindowResult/WindowResult.hpp"
#include "Model.hpp"
#include <vector>

class WindowModel : Model {
public:
  std::vector<Interval> windows, ref_intervals, query_intervals;
  ChrSizesVector chr_sizes;
  std::string method;
  Algorithm algorithm;

  WindowModel();
  WindowModel(std::vector<Interval> windows,
              std::vector<Interval> ref_intervals,
              std::vector<Interval> query_intervals, ChrSizesMap chr_sizes_map,
              Algorithm algorithm);

  std::vector<WindowResult> run();

  static std::vector<std::vector<Interval>>
  get_windows_intervals_naive(const std::vector<Interval> &windows,
                              const std::vector<Interval> &intervals);

  static std::vector<std::vector<Interval>>
  get_windows_intervals(const std::vector<Interval> &windows,
                        const std::vector<Interval> &intervals);

  std::vector<WindowResult> probs_by_window_single_chr(
      const std::vector<Interval> &windows,
      const std::vector<Interval> &windows_ref_intervals,
      const std::vector<Interval> &windows_query_intervals,
      const std::pair<std::string, long long> chr_size_entry);
};

#endif // WINDOWMODEL_H
