
#ifndef PROJECT_1_BRUTE_FORCE_NN_HPP
#define PROJECT_1_BRUTE_FORCE_NN_HPP

#include <map>
#include <vector>
#include <tuple>
#include "../util/dataset/dataset.hpp"

using std::multimap;
using std::vector;
using std::tuple;
using std::string;

// Computes NN with brute force
void bruteforce_nn(FlattenedCurve &query, vector<FlattenedCurve *> *data,
                    double(*distance_f)(FlattenedCurve&, FlattenedCurve&), std::tuple<double, string>* result);

#endif //PROJECT_1_BRUTE_FORCE_NN_HPP
