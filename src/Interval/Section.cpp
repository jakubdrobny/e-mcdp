#include "Section.hpp"

Section::Section(const std::string &chr_name, const long long &begin, const long long &end,
                 bool first_interval_intersected, bool last_interval_intersected)
    : chr_name(chr_name), begin(begin), end(end), first_interval_intersected(first_interval_intersected),
      last_interval_intersected(last_interval_intersected) {}

long long Section::length() { return this->end - this->begin; }
