
#include "lsh.hpp"
#include "../util/utilities.hpp"
#include <random>
#include <vector>
#include <list>
#include <map>
#include <algorithm>

using std::map;
using std::unordered_map;
using std::vector;
using std::list;
using std::string;
using std::get;

LSH::LSH(vector<FlattenedCurve *> *data, distance_f _distance, uint32_t num_ht, uint32_t num_hfs,
         uint32_t _k_nearest, double _radius)
: maps(new vector< hashtable *>()), distance(_distance), k_nearest_n(_k_nearest),
radius(_radius) {

    // Automatically compute the window size
    uint32_t window = estimate_window_size(data, this->distance);

    // Get the dimension of the data, by accessing a datapoint
    uint32_t dim = (*(*data)[0]).get_size();

    for(uint32_t i = 0; i < num_ht; i++) {
        /* Automatically scale the size of the Hash family's tables by computing:
            N / 2 ^ (log_10(N) - 1) */
        uint32_t size = data->size() / (uint32_t)pow(2, (log10((double)data->size()) - 1));
        this->maps->push_back(new hashtable(size, new amplified_hf(num_hfs, window, dim)));
    }

    // Load the data into the structure
    this->load(data);
}

LSH::~LSH() {
    for(auto ht: *this->maps)
        delete ht;
    this->maps->clear();
    delete this->maps;
}

/*
 *  Load data into the structure
 *  Process:
 *  1) Every amplified hash function is applied to a hashtable
 *  2) Run all points to all amplified hfs
 *  3) Result is L hash tables, each one including all points
 */
void LSH::load(vector<FlattenedCurve *> *data) {
    for(auto & map : *this->maps)
        for(auto pair: *data)
            map->insert(pair);
}

// Run K Nearest Neighbours
void LSH::knn(FlattenedCurve *query, multimap<double, string>& results) {

    // Save computed distances
    auto dist_cache = unordered_map<string, double>();
    // For every map
    for(auto & map : *this->maps) {
        // Find the bucket that the query is hashed into
        uint32_t id = map->hash(query);
        uint32_t index = id % map->get_tablesize();
        auto hashbucket = map->get_bucket(index);

        // For every point in the current bucket
        for(auto & it : *hashbucket) {

            uint32_t identity = get<0>(*it);
            FlattenedCurve *p = get<1>(*it);
            string label = p->get_id();

            /* Check for distances, iff item's key is equal to query's amplified hf value */
            if(id == identity) {
                // Compute the distance between the query and the current point
                double dist;
                if (dist_cache.find(label) == dist_cache.end()) {
                    dist = this->distance(*query, *p);
                    dist_cache.insert({label, dist});
                    // Store the appropriate info into the top N neighbours map
                    if (results.size() == this->k_nearest_n) {
                        if (dist < results.rbegin()->first)
                            results.erase(--results.end());
                        else continue;
                    }
                    results.insert({dist, label});
                }
            }
        }
    }
}

// Run Range-search
void LSH::range_search(FlattenedCurve *query, list<tuple<FlattenedCurve *, double>> &results) {

    // Save computed distances
    auto dist_cache = unordered_map<string, double>();
    // For every map
    for(auto & map : *this->maps) {
        // Find the bucket that the query is hashed into
        uint32_t id = map->hash(query);
        uint32_t index = id % map->get_tablesize();
        auto hashbucket = map->get_bucket(index);

        // For every point in the current bucket
        for(auto & it : *hashbucket) {

            uint32_t identity = get<0>(*it);
            FlattenedCurve *p = get<1>(*it);
            string label = p->get_id();

            /* Check for distances, iff items key is equal to query's amplified hf value */
            if(id == identity) {
                // Compute the distance between the query and the current point
                double dist;
                if (dist_cache.find(label) == dist_cache.end()) {
                    dist = this->distance(*query, *p);
                    dist_cache.insert({label, dist});
                    if (dist < this->radius)
                        results.emplace_back(p, dist);
                }
            }
        }
    }
}
