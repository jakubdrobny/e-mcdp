#include "WindowSectionSplitResult.hpp"
#include "../Interval/Interval.hpp"
#include <vector>

WindowSectionSplitResult::WindowSectionSplitResult() {}

WindowSectionSplitResult::WindowSectionSplitResult(std::vector<Section> sections, std::vector<Interval> spans)
    : sections(sections), spans(spans) {}

std::vector<Section> WindowSectionSplitResult::get_sections() { return sections; }
std::vector<Interval> WindowSectionSplitResult::get_spans() { return spans; }
