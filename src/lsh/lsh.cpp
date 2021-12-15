
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
         const std::string &_metric, uint32_t num_ht, uint32_t num_hfs,
         double _radius)
: maps(new vector< hashtable *>()),
  radius(_radius),
  grids(new std::vector<Grid *>()),
  raw_inputs(input.getData()),
  raw_queries(queries.getData()),
  L_flattened_inputs(new std::vector< std::vector< FlattenedCurve *> *>()),
  L_flattened_queries(new std::vector< std::vector< FlattenedCurve *> *>()) {

    set_metrics(_metric);

    // TODO: fine-tune this
    double grid_interval = estimate_grid_interval(input, queries);
    uint32_t max_curve_length = compute_max_curve_length(input, queries);

    if(_metric == "discr_frechet") {
        for (int i = 0; i < num_ht; ++i)
            grids->push_back(new Grid(grid_interval));

        max_curve_length *= 2;
        curves_preprocess(*input.getData(), max_curve_length, "input");
        curves_preprocess(*queries.getData(), max_curve_length, "query");
    }

    if(_metric == "cont_frechet") {

        grids->push_back(new Grid(grid_interval));

        // TODO: Μaybe compute it based in the avg point distances of curves
        double pruning_threshold = 10;
        curves_preprocess(*input.getData(), pruning_threshold, max_curve_length, "input");
        curves_preprocess(*queries.getData(), pruning_threshold, max_curve_length, "query");
    }
        // Automatically compute the window size
        uint32_t window = estimate_window_size(input.getData(), Metrics::euclidean);

        for (uint32_t i = 0; i < num_ht; i++) {
            /* Automatically scale the size of the Hash family's tables by computing:
                N / 2 ^ (log_10(N) - 1) */
            uint32_t size = raw_inputs->size() /
                            (uint32_t) pow(2, (log10((double) raw_inputs->size()) - 1));
            this->maps->push_back(new hashtable(size, new amplified_hf(num_hfs, window, max_curve_length)));
        }
        // Load the data into the structure
        this->load(L_flattened_inputs);
}

void LSH::set_metrics(const string &_metric) {

    if(_metric == "cont_frechet") {
        metric = (distance_f)Metrics::continuous_frechet_distance;
    }
    else if(_metric == "discr_frechet") {
        metric = (distance_f)Metrics::discrete_frechet_distance;
    }
    else {
        metric = (distance_f)Metrics::euclidean;
    }
}

void LSH::curves_preprocess(std::vector<Curve *> &curves, double pruning_threshold, uint32_t max_curve_length, const std::string& type) {
    int i = 0;
    for(auto & curve: curves) {
        for(auto &grid: *grids) {
            curve->filter(pruning_threshold);
            curve->erase_time_axis();
            grid->fit(*curve);
            curve->min_max_filter();
            //TODO: see if padding num is good
            curve->apply_padding(max_curve_length);
            auto flattened_curve = curve->flatten();

            if (type == "input") {
                if (L_flattened_inputs->size() < i + 1)
                    L_flattened_inputs->push_back(new std::vector<FlattenedCurve *>());
                (*L_flattened_inputs)[i]->push_back(flattened_curve);
            } else {
                if (L_flattened_queries->size() < i + 1)
                    L_flattened_queries->push_back(new std::vector<FlattenedCurve *>());
                (*L_flattened_queries)[i]->push_back(flattened_curve);
            }
        }
        i++;
    }
}

void LSH::curves_preprocess(std::vector<Curve *> &curves, uint32_t max_curve_length, const std::string& type) {

    int i = 0;
    for(auto &curve: curves) {
        for(auto &grid: *grids) {

            auto _curve = *curve;
            grid->fit(_curve);
            grid->remove_consecutive_duplicates(_curve);
            auto flattened_curve = _curve.flatten();
            flattened_curve->apply_padding(max_curve_length);

            if(type == "input") {
                if(L_flattened_inputs->size() < i+1)
                    L_flattened_inputs->push_back(new std::vector<FlattenedCurve *>());
                (*L_flattened_inputs)[i]->push_back(flattened_curve);
            }
            else {
                if(L_flattened_queries->size() < i+1)
                    L_flattened_queries->push_back(new std::vector<FlattenedCurve *>());
                (*L_flattened_queries)[i]->push_back(flattened_curve);
            }
        }
        i++;
    }
}

void LSH::delete_flattened_inputs() {
    for (auto &curve: *L_flattened_inputs) {
        for (int i = 0; i < curve->size(); i++) {
            delete curve->data()[i];
        }
        delete curve;
    }
    delete L_flattened_inputs;

    for (auto &curve: *L_flattened_queries) {
        for (int i = 0; i < curve->size(); i++) {
            delete curve->data()[i];
        }
        delete curve;
    }
    delete L_flattened_queries;
}

LSH::~LSH() {

    for(auto ht: *this->maps)
        delete ht;
    this->maps->clear();
    delete this->maps;

    for(auto &grid: *grids)
        delete grid;
    delete grids;

    delete_flattened_inputs();
}

/*
 *  Load data into the structure
 *  Process:
 *  1) Every amplified hash function is applied to a hashtable
 *  2) Run all points to all amplified hfs
 *  3) Result is L hash tables, each one including all points
 */
void LSH::load(vector< vector<FlattenedCurve *> *> *data) {
    int i = 0;
    for(auto & map : *this->maps) {
        for (auto curve_group: *data) {
            map->insert((*curve_group)[i]);
        }
        i++;
    }
}

// Run K Nearest Neighbours
void LSH::nn(vector<FlattenedCurve *> &query_family, std::tuple<double, string> &result) {

    // Save computed distances
    auto dist_cache = unordered_map<string, double>();
    // For every map
    int i = 0;
    for(auto & map : *this->maps) {
        // Find the bucket that the query is hashed into
        uint32_t id = map->hash(query_family[i]);
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
                    dist = metric(*query_family[i], *p);
                    dist_cache.insert({label, dist});

                    if(dist < get<0>(result)) {
                        std::get<0>(result) = dist;
                        std::get<1>(result) = label;
                    }
                }
            }
        }
        i++;
    }
}

void LSH::nearest_neighbor(const std::string &out_path) {

    double avg_lsh_time_taken = 0;
    double abg_brutef_time_taken = 0;
    double maf = 0;

    int i = 0;
    for(auto query_family: *L_flattened_queries) {

        string label = (*query_family)[0]->get_id();

        // Benchmark LSH K-NN
        auto start = GET_CURR_TIME();
        std::tuple<double, string> top_lsh = {std::numeric_limits<double>::max(), "-"};
        nn(*query_family, top_lsh);
        auto end = GET_CURR_TIME();
        double lsh_time_taken = GET_DURATION(start, end);

        // Benchmark Brute-force K-NN
        start = GET_CURR_TIME();
        std::tuple<double, string> top_brutef = {std::numeric_limits<double>::max(), "-"};
        bruteforce_nn(*(*raw_queries)[i], raw_inputs, Metrics::euclidean, &top_brutef);
        end = GET_CURR_TIME();
        auto brutef_time_taken = GET_DURATION(start, end);

        avg_lsh_time_taken += (lsh_time_taken / raw_queries->size());
        abg_brutef_time_taken += (brutef_time_taken / raw_queries->size());
        double candidate_maf = std::get<0>(top_lsh) / std::get<0>(top_brutef);
        if(maf < candidate_maf)
            maf = candidate_maf;

        write_data_to_out_file(label, top_lsh, top_brutef, out_path);

        i++;
    }
    write_data_to_out_file(avg_lsh_time_taken, abg_brutef_time_taken, maf, out_path);
}

// Run Range-search
/*void LSH::range_search(FlattenedCurve *query, list<tuple<FlattenedCurve *, double>> &results) {

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

            *//* Check for distances, iff items key is equal to query's amplified hf value *//*
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
}*/
