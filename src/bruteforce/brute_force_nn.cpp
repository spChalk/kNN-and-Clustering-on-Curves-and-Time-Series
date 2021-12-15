
#include <iostream>
#include "brute_force_nn.hpp"

using std::get;
using std::multimap;

void bruteforce_nn(Curve &query, vector<Curve *> *data,
                                    double(*distance_f)(Curve&, Curve&), std::tuple<double, string>* result) {
    for(auto curve: *data) {
        string label = curve->get_id();
        // Compute the distance between the query and the current point
        double dist = distance_f(query, *curve);
        if(dist < get<0>(*result)) {
            std::get<0>(*result) = dist;
            std::get<1>(*result) = label;
        }
    }
}
