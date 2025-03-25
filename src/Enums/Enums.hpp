#ifndef ENUM_H
#define ENUM_H

#include <map>
#include <string>

enum class Algorithm { NAIVE, FAST };

template <typename T> extern bool validate_enum(const std::map<std::string, T> &stringToEnum, const std::string &str) {
  return stringToEnum.count(str);
}

extern const std::map<std::string, Algorithm> algorithmToEnum;
extern const std::map<Algorithm, std::string> algorithmToString;

#endif // ENUM_H
