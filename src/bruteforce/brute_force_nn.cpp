
#include <iostream>
#include "brute_force_nn.hpp"

using std::get;
using std::multimap;

void bruteforce_nn(FlattenedCurve &query, vector<FlattenedCurve *> *data,
                                    double(*distance_f)(FlattenedCurve&, FlattenedCurve&), std::tuple<double, string>* result) {
    for(auto pair: *data) {
        string label = pair->get_id();
        // Compute the distance between the query and the current point
        double dist = distance_f(query, *pair);
        if(dist < get<0>(*result)) {
            std::get<0>(*result) = dist;
            std::get<1>(*result) = label;
        }
    }
}