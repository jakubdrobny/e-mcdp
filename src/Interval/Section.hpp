#ifndef SECTION_H
#define SECTION_H

#include "../Results/SectionProbs.hpp"
#include "Interval.hpp"
#include <string>

// half open intervals [b, e), but with extra info (whether first interval is from the previous section and wheter last
// interval is going into the next section)
class Section : public Interval {
public:
  Section(const std::string &chr_name, const long long &begin, const long long &end, bool first_interval_intersected,
          bool last_interval_intersected);
  Section(const std::string &chr_name, const long long &begin, const long long &end, bool first_interval_intersected,
          bool last_interval_intersected, const std::vector<Interval> &intervals);

  bool get_first_interval_intersected() const;
  bool get_last_interval_intersected() const;
  std::vector<Interval> get_intervals() const;
  void set_intervals(const std::vector<Interval> &new_intervals);
  SectionProbs get_probs() const;
  void set_probs(const SectionProbs &new_probs);
  long long get_overlap_count() const;
  void set_overlap_count(long long new_overlap_count);

private:
  bool first_interval_intersected, last_interval_intersected;
  std::vector<Interval> intervals;
  SectionProbs probs;
  long long overlap_count;
};

#endif // INTERVAL_H
