
#ifndef PROJECT_1_BRUTE_FORCE_KNN_HPP
#define PROJECT_1_BRUTE_FORCE_KNN_HPP

#include <map>
#include <vector>
#include <tuple>
#include "../util/dataset/dataset.hpp"

using std::multimap;
using std::vector;
using std::tuple;
using std::string;

// Computes K-NN with brute force
void bruteforce_knn(uint32_t k, vector<double> &query, vector<Point *> *data,
                                    double(*distance_f)(Point&, Point&), multimap<double, string>*);

#endif //PROJECT_1_BRUTE_FORCE_KNN_HPP
