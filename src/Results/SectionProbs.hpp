#ifndef SECTIONPROBS_H
#define SECTIONPROBS_H

#include "WindowResult.hpp"

class SectionProbs {
public:
  SectionProbs();
  SectionProbs(MultiProbs normal, MultiProbs except_first, MultiProbs except_last, MultiProbs except_first_and_last);

  MultiProbs get_normal() const;
  MultiProbs get_except_first() const;
  MultiProbs get_except_last() const;
  MultiProbs get_except_first_and_last() const;

  SectionProbs operator*(const SectionProbs &other) const;

private:
  MultiProbs normal, except_first, except_last, except_first_and_last;
};

#endif // SECTIONPROBS_H
