#ifndef WINDOWSECTIONSPLITRESULT_H
#define WINDOWSECTIONSPLITRESULT_H

#include "../Interval/Interval.hpp"
#include <vector>

class WindowSectionSplitResult {
private:
  std::vector<Interval> sections, spans;

public:
  WindowSectionSplitResult();
  WindowSectionSplitResult(std::vector<Interval> sections,
                           std::vector<Interval> spans);

  std::vector<Interval> get_sections();
  std::vector<Interval> get_spans();
};

#endif // WINDOWSECTIONSPLITRESULT_H
