#include "Enums.hpp"

const std::map<std::string, Algorithm> algorithmToEnum = {{"naive", Algorithm::NAIVE},
                                                          {"slow", Algorithm::SLOW},
                                                          {"slow_bad", Algorithm::SLOW_BAD},
                                                          {"fast", Algorithm::FAST},
                                                          {"fast", Algorithm::FAST_BAD}};
const std::map<std::string, Statistic> statisticToEnum = {{"overlaps", Statistic::OVERLAPS},
                                                          {"bases", Statistic::BASES}};

const std::map<Statistic, std::string> statisticToString = {{Statistic::OVERLAPS, "overlaps"},
                                                            {Statistic::BASES, "bases"}};
const std::map<Algorithm, std::string> algorithmToString = {{Algorithm::NAIVE, "naive"},
                                                            {Algorithm::SLOW, "slow"},
                                                            {Algorithm::SLOW_BAD, "slow_bad"},
                                                            {Algorithm::FAST, "fast"},
                                                            {Algorithm::FAST_BAD, "fast_bad"}};
