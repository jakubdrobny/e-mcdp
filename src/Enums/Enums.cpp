#include "Enums.hpp"

const std::map<std::string, Algorithm> algorithmToEnum = {
    {"naive", Algorithm::NAIVE}, {"fast", Algorithm::FAST}};
const std::map<Algorithm, std::string> algorithmToString = {
    {Algorithm::NAIVE, "naive"}, {Algorithm::FAST, "fast"}};
