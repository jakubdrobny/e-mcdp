#ifndef WINDOWMODEL_H
#define WINDOWMODEL_H

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

  WindowModel();
  WindowModel(std::vector<Interval> windows,
              std::vector<Interval> ref_intervals,
              std::vector<Interval> query_intervals, ChrSizesMap chr_sizes_map);

  std::vector<WindowResult> run();

  static std::vector<std::vector<Interval>>
  get_windows_intervals(const std::vector<Interval> &windows,
                        const std::vector<Interval> &intervals);
};

#endif // WINDOWMODEL_H
