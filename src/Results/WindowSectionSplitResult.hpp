#ifndef WINDOWSECTIONSPLITRESULT_H
#define WINDOWSECTIONSPLITRESULT_H

#include "../Interval/Interval.hpp"
#include "../Interval/Section.hpp"
#include <vector>

class WindowSectionSplitResult {
private:
  std::vector<Section> sections;
  std::vector<Interval> spans;

public:
  WindowSectionSplitResult();
  WindowSectionSplitResult(std::vector<Section> sections, std::vector<Interval> spans);

  std::vector<Section> get_sections();
  std::vector<Interval> get_spans();
};

#endif // WINDOWSECTIONSPLITRESULT_H
