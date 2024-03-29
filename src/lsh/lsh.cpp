
#include "lsh.hpp"
#include "../bruteforce/brute_force_nn.hpp"
#include "../util/files/file_reader.hpp"
#include <random>
#include <unordered_map>
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

#define PRUNING_THRESHOLD (10)
#define WINDOW_SIZE (1000)

LSH::LSH(Dataset &input, const enum metrics &_metric, uint32_t num_ht, uint32_t num_hfs, double _radius):
        maps(new vector< hashtable *>()),
        metric_id(_metric),
        radius(_radius),
        grids(new std::vector<Grid *>()),
        raw_inputs(input.getData()),
        raw_queries(nullptr),
        padding_len(compute_max_curve_length(input.getData())),
        L_flattened_inputs(new vector_of_flattened_curves()),
        L_flattened_queries(new vector_of_flattened_curves()),
        label_to_curve(new std::unordered_map<std::string, Curve *>()),
        label_to_index_in_fl_queries(new std::unordered_map<std::string, uint32_t>()) {

    if(metric_id == DISCRETE_FRECHET) padding_len *= 2;

    for (uint32_t i = 0; i < num_ht; i++) {
        /* Automatically scale the size of the Hash family's tables by computing:
            N / 2 ^ (log_10(N) - 1) */
        uint32_t size = raw_inputs->size() /
                        (uint32_t) pow(2, (log10((double) raw_inputs->size()) - 1));
        this->maps->push_back(new hashtable(size, new amplified_hf(num_hfs, WINDOW_SIZE, padding_len)));
    }

    set_metrics_and_preprocess();

    // Load the data into the structure
    this->load();
}

LSH::LSH(vector<FlattenedCurve *> *data, const enum metrics &_metric, uint32_t num_ht, uint32_t num_hfs,
         uint32_t _k_nearest, double _radius)
: maps(new vector< hashtable *>()), metric_id(_metric), radius(_radius) {

    // Automatically compute the window size
    uint32_t window = WINDOW_SIZE;

    // Get the dimension of the data, by accessing a datapoint
    uint32_t dim = data->at(0)->get_size();

    for(uint32_t i = 0; i < num_ht; i++) {
        /* Automatically scale the size of the Hash family's tables by computing:
            N / 2 ^ (log_10(N) - 1) */
        uint32_t size = data->size() / (uint32_t)pow(2, (log10((double)data->size()) - 1));
        this->maps->push_back(new hashtable(size, new amplified_hf(num_hfs, window, dim)));
    }

    // Load the data into the structure
    for(auto & map : *this->maps)
        for(auto pair: *data)
            map->insert(pair);
}

// Set the appropriate metrics and preprocess the input data
void LSH::set_metrics_and_preprocess() {

    double grid_interval = estimate_grid_interval(raw_inputs);

    if(metric_id == EUCLIDEAN)
        metric = Metrics::Euclidean::distance;
    else {
        for (uint32_t i = 0; i < maps->size(); ++i)
            grids->push_back(new Grid(grid_interval));

        metric = metric_id == CONTINUOUS_FRECHET ? (distance_f)Metrics::Continuous_Frechet::distance :
                 Metrics::Discrete_Frechet::distance;
    }

    curves_preprocess(*raw_inputs, "input");
}

void LSH::curves_preprocess(std::vector<Curve *> &curves, const std::string& type) {
    for(auto & curve: curves)
        curve_preprocess(*curve, type);
}

void LSH::curve_preprocess(Curve &c, const string &type) {

    label_to_curve->insert({c.get_id(), &c});

    if(type == "input")
        L_flattened_inputs->push_back(new flattened_curves());
    else {
        L_flattened_queries->push_back(new flattened_curves());
        // Keep the query's index in order to remember it later, if needed.
        label_to_index_in_fl_queries->insert({c.get_id(), L_flattened_queries->size() - 1});
    }

    for(uint32_t j = 0; j < maps->size(); ++j) {

        auto flattened_curve = metric_id == EUCLIDEAN ? euclidean_preprocess(c) :
                               metric_id == CONTINUOUS_FRECHET ? cont_frechet_preprocess(c, j) :
                               discr_frechet_preprocess(c, j);

        if(type == "input")
            L_flattened_inputs->back()->push_back(flattened_curve);
        else
            L_flattened_queries->back()->push_back(flattened_curve);
    }
}

FlattenedCurve *LSH::euclidean_preprocess(Curve &curve) {
    auto _curve = curve;
    _curve.erase_time_axis();
    return _curve.flatten();
}

FlattenedCurve *LSH::euclidean_preprocess(FlattenedCurve &curve) {
    return new FlattenedCurve(curve);
}

FlattenedCurve *LSH::cont_frechet_preprocess(Curve &curve, uint32_t index) {
    auto _curve = curve;
    _curve.filter(PRUNING_THRESHOLD);
    _curve.erase_time_axis();
    (*grids)[index]->fit(_curve);
    _curve.min_max_filter();
    _curve.apply_padding(padding_len);
    return _curve.flatten();
}

FlattenedCurve *LSH::discr_frechet_preprocess(Curve &curve, uint32_t index) {
    auto _curve = Curve(curve);
    (*grids)[index]->fit(_curve);
    (*grids)[index]->remove_consecutive_duplicates(_curve);
    auto *flattened_curve = _curve.flatten();
    flattened_curve->apply_padding(padding_len);
    return flattened_curve;
}

void LSH::delete_flattened_inputs() {
    if (L_flattened_inputs) {
        for (auto &curve: *L_flattened_inputs) {
            for (uint32_t i = 0; i < curve->size(); i++) {
                delete curve->data()[i];
            }
            delete curve;
        }
        delete L_flattened_inputs;
    }

    if (L_flattened_queries) {
        for (auto &curve: *L_flattened_queries) {
            for (uint32_t i = 0; i < curve->size(); i++) {
                delete curve->data()[i];
            }
            delete curve;
        }
        delete L_flattened_queries;
    }
}

LSH::~LSH() {

    for(auto ht: *this->maps)
        delete ht;
    this->maps->clear();
    delete this->maps;

    if (grids) {
        for(auto &grid: *grids)
            delete grid;
        delete grids;
    }

    delete_flattened_inputs();

    if (label_to_curve) {
        label_to_curve->clear();
        delete label_to_curve;
    }

    if (label_to_index_in_fl_queries) {
        label_to_index_in_fl_queries->clear();
        delete label_to_index_in_fl_queries;
    }
}

/*
 *  Load data into the structure
 *  Process:
 *  1) Every amplified hash function is applied to a hashtable
 *  2) Run all points to all amplified hfs
 *  3) Result is L hash tables, each one including all points
 */
void LSH::load() {
    int i = 0;
    for(auto & map : *this->maps) {
        for (auto curve_group: *L_flattened_inputs)
            map->insert((*curve_group)[i]);
        i++;
    }
}

// Run K Nearest Neighbours
void LSH::nn(flattened_curves &query_family, std::tuple<double, string> &result) {

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

                    if(metric_id != EUCLIDEAN) {
                        Curve *raw_query_curve = label_to_curve->find((*query_family[i]).get_id())->second;
                        Curve *raw_current_curve = label_to_curve->find(label)->second;
                        dist = metric(*raw_query_curve, *raw_current_curve);
                    } else
                        dist = Metrics::Euclidean::distance(*query_family[i], *p);

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

// Receives a label and returns the corresponding query's gridded family.
flattened_curves *LSH::get_flattened_family(std::string &label) {
    uint32_t index = label_to_index_in_fl_queries->find(label)->second;
    assert(index <= L_flattened_queries->size());
    return (*L_flattened_queries)[index];
}

void LSH::nearest_neighbor(Curve *query, std::tuple<double, string> &result) {

    // Preprocess query, if needed
    curve_preprocess(*query, "query");
    auto label = query->get_id();
    // Get the query's flattened family from L_flattened_queries
    auto query_family = get_flattened_family(label);

    // Run NN
    nn(*query_family, result);
}


// Original Range Search from HW1
void LSH::range_search(FlattenedCurve *query, std::list<std::tuple<FlattenedCurve *, double>> & results)
{
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
                    dist = Metrics::Euclidean::distance(*query, *p);
                    dist_cache.insert({label, dist});
                    if (dist < this->radius)
                        results.emplace_back(p, dist);
                }
            }
        }
    }
}


void LSH::range_search(Curve *query, list<tuple<Curve *, double>> &results)
{
    curve_preprocess(*query, "query");
    auto label = query->get_id();
    // Get the query's flattened family from L_flattened_queries
    auto *query_family = get_flattened_family(label);

    // Save computed distances
    auto dist_cache = unordered_map<string, double>();
    // For every map
    int i = 0;
    for(auto & map : *this->maps) {
        // Find the bucket that the query is hashed into
        uint32_t id = map->hash(query_family->at(i));
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

                    Curve *raw_current_curve = label_to_curve->find(label)->second;

                    if(metric_id != EUCLIDEAN) {
                        Curve *raw_query_curve = label_to_curve->find((*query_family->at(i)).get_id())->second;
                        dist = metric(*raw_query_curve, *raw_current_curve);

                    } else
                        dist = Metrics::Euclidean::distance(*query_family->at(i), *p);

                    dist_cache.insert({label, dist});

                    if (dist < this->radius)
                        results.emplace_back(raw_current_curve, dist);
                }
            }
        }
        ++i;
    }
}