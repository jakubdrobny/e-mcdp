#ifndef SECTION_H
#define SECTION_H

#include "../Results/SectionProbs.hpp"
#include "Interval.hpp"
#include <string>

// half open intervals [b, e), but with extra info (whether first interval is from the previous section and wheter last
// interval is going into the next section)
class Section : public Interval {
public:
  Section(const std::string &chr_name, const long long &begin, const long long &end,
          bool first_ref_interval_intersected, bool last_ref_interval_intersected,
          bool first_query_interval_intersected, bool last_query_interval_intersected);
  Section(const std::string &chr_name, const long long &begin, const long long &end, bool first_interval_intersected,
          bool last_ref_interval_intersected, bool first_query_interval_intersected,
          bool last_query_interval_intersected, const std::vector<Interval> &ref_intervals,
          const std::vector<Interval> &query_intervals);

  bool get_first_ref_interval_intersected() const;
  bool get_last_ref_interval_intersected() const;
  bool get_first_query_interval_intersected() const;
  bool get_last_query_interval_intersected() const;
  std::vector<Interval> get_ref_intervals() const;
  void set_ref_intervals(const std::vector<Interval> &new_ref_intervals);
  std::vector<Interval> get_query_intervals() const;
  void set_query_intervals(const std::vector<Interval> &new_query_intervals);
  SectionProbs get_probs() const;
  void set_probs(const SectionProbs &new_probs);
  long long get_overlap_count() const;
  void set_overlap_count(long long new_overlap_count);

private:
  bool first_ref_interval_intersected, last_ref_interval_intersected, first_query_interval_intersected,
      last_query_interval_intersected;
  std::vector<Interval> ref_intervals, query_intervals;
  SectionProbs probs;
  long long overlap_count;
};

std::ostream &operator<<(std::ostream &os, const Section &section);

#endif // INTERVAL_H
