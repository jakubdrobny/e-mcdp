#include "Enums.hpp"

const std::map<std::string, Algorithm> algorithmToEnum = {
    {"naive", Algorithm::NAIVE}, {"fast", Algorithm::FAST}, {"test", Algorithm::TEST}};
const std::map<std::string, Statistic> statisticToEnum = {{"overlaps", Statistic::OVERLAPS},
                                                          {"bases", Statistic::BASES}};

const std::map<Statistic, std::string> statisticToString = {{Statistic::OVERLAPS, "overlaps"},
                                                            {Statistic::BASES, "bases"}};
const std::map<Algorithm, std::string> algorithmToString = {
    {Algorithm::NAIVE, "naive"}, {Algorithm::FAST, "fast"}, {Algorithm::TEST, "test"}};
