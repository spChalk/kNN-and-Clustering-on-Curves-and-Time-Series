
#include "lsh.hpp"
#include "../util/utilities.hpp"
#include "../util/metrics/metrics.hpp"
#include "../bruteforce/brute_force_nn.hpp"
#include "../util/files/file_reader.hpp"
#include <random>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <chrono>

using std::map;
using std::unordered_map;
using std::vector;
using std::list;
using std::string;
using std::get;

#define GET_DURATION(START, END) (std::chrono::duration_cast<std::chrono::nanoseconds>((END) - (START)).count() * 1e-9)
#define GET_CURR_TIME() (std::chrono::high_resolution_clock::now())

LSH::LSH(Dataset &input, Dataset &queries,
         distance_f _distance, uint32_t num_ht, uint32_t num_hfs,
         double _radius)
: maps(new vector< hashtable *>()),
  distance(_distance),
  radius(_radius),
  grid(nullptr),
  flattened_input_data(nullptr),
  flattened_query_data(nullptr) {

    // TODO: Μaybe compute it based in the avg point distances of curves
    double pruning_threshold = 10;

    // TODO: fine-tune this
    double grid_interval = estimate_grid_interval(input, queries);
    uint32_t max_curve_length = compute_max_curve_length(input, queries);

    grid = new Grid(grid_interval, (*input.getData())[0]->get_data_dimensions());

    curves_preprocess(*input.getData(), pruning_threshold, max_curve_length);
    curves_preprocess(*queries.getData(), pruning_threshold, max_curve_length);

    flattened_input_data = input.flatten_data();
    flattened_query_data = queries.flatten_data();

    // Automatically compute the window size
    uint32_t window = estimate_window_size(flattened_input_data, this->distance);

    // Get the dimension of the data, by accessing a datapoint
    uint32_t dim = (*(*flattened_input_data)[0]).get_size();

    for(uint32_t i = 0; i < num_ht; i++) {
        /* Automatically scale the size of the Hash family's tables by computing:
            N / 2 ^ (log_10(N) - 1) */
        uint32_t size = flattened_input_data->size() / (uint32_t)pow(2, (log10((double)flattened_input_data->size()) - 1));
        this->maps->push_back(new hashtable(size, new amplified_hf(num_hfs, window, dim)));
    }

    // Load the data into the structure
    this->load(flattened_input_data);
}

void LSH::curves_preprocess(std::vector<Curve *> &curves, double pruning_threshold, uint32_t max_curve_length) {
    for(auto & curve: curves) {
        curve->filter(pruning_threshold);
        curve->erase_time_axis();
        grid->fit(*curve);
        curve->min_max_filter();
        //TODO: see if padding num is good
        curve->apply_padding(max_curve_length);
    }
}

LSH::~LSH() {
    for(auto ht: *this->maps)
        delete ht;
    this->maps->clear();
    delete this->maps;

    for (auto &d: *flattened_input_data)
        delete d;
    delete flattened_input_data;
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
void LSH::nn(FlattenedCurve &query, std::tuple<double, string> &result) {

    // Save computed distances
    auto dist_cache = unordered_map<string, double>();
    // For every map
    for(auto & map : *this->maps) {
        // Find the bucket that the query is hashed into
        uint32_t id = map->hash(&query);
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
                    dist = this->distance(query, *p);
                    dist_cache.insert({label, dist});

                    if(dist < get<0>(result)) {
                        std::get<0>(result) = dist;
                        std::get<1>(result) = label;
                    }
                }
            }
        }
    }
}

void LSH::nearest_neighbor(const std::string &out_path) {

    double avg_lsh_time_taken = 0;
    double abg_brutef_time_taken = 0;
    double maf = 0;

    for(auto curve: *flattened_query_data) {

        string label = curve->get_id();

        // Benchmark LSH K-NN
        auto start = GET_CURR_TIME();
        std::tuple<double, string> top_lsh = {std::numeric_limits<double>::max(), "-"};
        nn(*curve, top_lsh);
        auto end = GET_CURR_TIME();
        double lsh_time_taken = GET_DURATION(start, end);

        // Benchmark Brute-force K-NN
        start = GET_CURR_TIME();
        std::tuple<double, string> top_brutef = {std::numeric_limits<double>::max(), "-"};
        bruteforce_nn(*curve, flattened_input_data, Metrics::euclidean, &top_brutef);
        end = GET_CURR_TIME();
        auto brutef_time_taken = GET_DURATION(start, end);

        avg_lsh_time_taken += (lsh_time_taken / flattened_query_data->size());
        abg_brutef_time_taken += (brutef_time_taken / flattened_query_data->size());
        double candidate_maf = std::get<0>(top_lsh) / std::get<0>(top_brutef);
        if(maf < candidate_maf)
            maf = candidate_maf;

        write_data_to_out_file(label, top_lsh, top_brutef, out_path);
    }
    write_data_to_out_file(avg_lsh_time_taken, abg_brutef_time_taken, maf, out_path);
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
