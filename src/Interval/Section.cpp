#include "Section.hpp"
#include <vector>

Section::Section(const std::string &chr_name, const long long &begin, const long long &end,
                 bool first_interval_intersected, bool last_interval_intersected)
    : Interval(chr_name, begin, end), first_interval_intersected(first_interval_intersected),
      last_interval_intersected(last_interval_intersected) {}

Section::Section(const std::string &chr_name, const long long &begin, const long long &end,
                 bool first_interval_intersected, bool last_interval_intersected,
                 const std::vector<Interval> &intervals)
    : Interval(chr_name, begin, end), first_interval_intersected(first_interval_intersected),
      last_interval_intersected(last_interval_intersected), intervals(intervals) {}

bool Section::get_first_interval_intersected() const { return this->first_interval_intersected; }

bool Section::get_last_interval_intersected() const { return this->last_interval_intersected; }

std::vector<Interval> Section::get_intervals() const { return this->intervals; }

void Section::set_intervals(const std::vector<Interval> &new_intervals) { this->intervals = new_intervals; }

SectionProbs Section::get_probs() const { return this->probs; }

void Section::set_probs(const SectionProbs &new_probs) { this->probs = new_probs; }

long long Section::get_overlap_count() const { return this->overlap_count; }

void Section::set_overlap_count(long long new_overlap_count) { this->overlap_count = new_overlap_count; }
