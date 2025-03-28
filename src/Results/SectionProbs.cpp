#include "SectionProbs.hpp"

SectionProbs::SectionProbs() {}
SectionProbs::SectionProbs(MultiProbs normal, MultiProbs except_first, MultiProbs except_last,
                           MultiProbs except_first_and_last)
    : normal(normal), except_first(except_first), except_last(except_last),
      except_first_and_last(except_first_and_last) {}

MultiProbs SectionProbs::get_normal() const { return this->normal; }
MultiProbs SectionProbs::get_except_first() const { return this->except_first; }
MultiProbs SectionProbs::get_except_last() const { return this->except_last; }
MultiProbs SectionProbs::get_except_first_and_last() const { return this->except_first_and_last; }
