
#include <iostream>
#include "brute_force_knn.hpp"

using std::get;
using std::multimap;

void bruteforce_knn(uint32_t k, vector<double> &query, vector<FlattenedCurve *> *data,
                                    double(*distance_f)(FlattenedCurve&, FlattenedCurve&), multimap<double, string>* results) {
    if(k == 0) return;
    for(auto pair: *data) {
        
        string label = pair->get_id();

        // Compute the distance between the query and the current point
        std::string s;
        auto q = FlattenedCurve(s, query);
        double dist = distance_f(q, *pair);

        // Store the appropriate info into the top N neighbours map
        if(results->size() == k) {
            if(dist < results->rbegin()->first)
                results->erase(--results->end());
            else continue;
        }
        results->insert({dist, label});
    }
}