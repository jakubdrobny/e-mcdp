#ifndef WINDOWMODEL_H
#define WINDOWMODEL_H

#include "../Enums/Enums.hpp"
#include "../Helpers/Helpers.hpp"
#include "../Interval/Interval.hpp"
#include "../Results/SectionProbs.hpp"
#include "../Results/WindowResult.hpp"
#include "Model.hpp"
#include <vector>

class WindowModel : Model {
public:
  std::vector<Interval> windows, ref_intervals, query_intervals;
  ChrSizesVector chr_sizes;
  std::string method;
  Algorithm algorithm;

  WindowModel();
  WindowModel(std::vector<Interval> windows, std::vector<Interval> ref_intervals, std::vector<Interval> query_intervals,
              ChrSizesMap chr_sizes_map, Algorithm algorithm);

  std::vector<WindowResult> run();

  static std::vector<std::vector<Interval>> get_windows_intervals_naive(const std::vector<Interval> &windows,
                                                                        const std::vector<Interval> &intervals);

  template <typename WindowType>
  static std::vector<std::vector<Interval>> get_windows_intervals(const std::vector<WindowType> &windows,
                                                                  const std::vector<Interval> &intervals);

  std::vector<WindowResult> probs_by_window_single_chr_naive(const std::vector<Interval> &windows,
                                                             const std::vector<Interval> &windows_ref_intervals,
                                                             const std::vector<Interval> &windows_query_intervals,
                                                             const std::pair<std::string, long long> chr_size_entry);

  std::vector<WindowResult> probs_by_window_single_chr_smarter(const std::vector<Interval> &windows,
                                                               const std::vector<Interval> &windows_ref_intervals,
                                                               const std::vector<Interval> &windows_query_intervals,
                                                               const std::pair<std::string, long long> chr_size_entry,
                                                               bool use_segtree = true);

  std::vector<WindowResult> probs_by_window_single_chr_smarter_new(
      const std::vector<Interval> &windows, const std::vector<Interval> &ref_intervals,
      const std::vector<Interval> &query_intervals, const std::pair<std::string, long long> chr_size_entry,
      bool use_segtree = true);

private:
  SectionProbs eval_probs_single_section(const Section &section, const MarkovChain &markov_chain);

  SectionProbs eval_probs_single_section_new(const Section &section, const MarkovChain &markov_chain);

  void correct_ends(Section &section, const MarkovChain &markov_chain);
};

#endif // WINDOWMODEL_H
