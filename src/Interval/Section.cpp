#include "Section.hpp"
#include <vector>

Section::Section(const std::string &chr_name, const long long &begin, const long long &end,
                 bool first_ref_interval_intersected, bool last_ref_interval_intersected,
                 bool first_query_interval_intersected, bool last_query_interval_intersected)
    : Interval(chr_name, begin, end), first_ref_interval_intersected(first_ref_interval_intersected),
      last_ref_interval_intersected(last_ref_interval_intersected),
      first_query_interval_intersected(first_query_interval_intersected),
      last_query_interval_intersected(last_query_interval_intersected) {}

Section::Section(const std::string &chr_name, const long long &begin, const long long &end,
                 bool first_ref_interval_intersected, bool last_ref_interval_intersected,
                 bool first_query_interval_intersected, bool last_query_interval_intersected,
                 const std::vector<Interval> &ref_intervals, const std::vector<Interval> &query_intervals)
    : Interval(chr_name, begin, end), first_ref_interval_intersected(first_ref_interval_intersected),
      last_ref_interval_intersected(last_ref_interval_intersected),
      first_query_interval_intersected(first_query_interval_intersected),
      last_query_interval_intersected(last_query_interval_intersected), ref_intervals(ref_intervals),
      query_intervals(query_intervals) {}

bool Section::get_first_ref_interval_intersected() const { return this->first_ref_interval_intersected; }

bool Section::get_last_ref_interval_intersected() const { return this->last_ref_interval_intersected; }

bool Section::get_first_query_interval_intersected() const { return this->first_query_interval_intersected; }

bool Section::get_last_query_interval_intersected() const { return this->last_query_interval_intersected; }

std::vector<Interval> Section::get_ref_intervals() const { return this->ref_intervals; }

void Section::set_ref_intervals(const std::vector<Interval> &new_ref_intervals) {
  this->ref_intervals = new_ref_intervals;
}

std::vector<Interval> Section::get_query_intervals() const { return this->query_intervals; }

void Section::set_query_intervals(const std::vector<Interval> &new_query_intervals) {
  this->query_intervals = new_query_intervals;
}

SectionProbs Section::get_probs() const { return this->probs; }

void Section::set_probs(const SectionProbs &new_probs) { this->probs = new_probs; }

long long Section::get_overlap_count() const { return this->overlap_count; }

void Section::set_overlap_count(long long new_overlap_count) { this->overlap_count = new_overlap_count; }
