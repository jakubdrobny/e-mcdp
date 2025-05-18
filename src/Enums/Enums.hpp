#ifndef ENUM_H
#define ENUM_H

#include <map>
#include <string>

enum class Algorithm { NAIVE, SLOW_BAD, SLOW, FAST_BAD, FAST };
enum class Statistic { OVERLAPS, BASES };
enum class Significance { ENRICHMENT, DEPLETION, COMBINED };

template <typename T> extern bool validate_enum(const std::map<std::string, T> &stringToEnum, const std::string &str) {
  return stringToEnum.count(str);
}

extern const std::map<std::string, Algorithm> algorithmToEnum;
extern const std::map<std::string, Statistic> statisticToEnum;
extern const std::map<std::string, Significance> significanceToEnum;

extern const std::map<Algorithm, std::string> algorithmToString;
extern const std::map<Statistic, std::string> statisticToString;
extern const std::map<Significance, std::string> significanceToString;

#endif // ENUM_H
