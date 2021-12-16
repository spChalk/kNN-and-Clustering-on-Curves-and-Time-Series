
#include "lsh.hpp"
#include "../bruteforce/brute_force_nn.hpp"
#include "../util/files/file_reader.hpp"
#include <random>
#include <unordered_map>
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

// TODO: Tune PRUNING_THRESHOLD (tuning threshold is used in time-series' dim. reduction)
// TODO: Fine-tune grid interval (δ)
// TODO: See if padding num is good (currently its max len of input data)
// TODO: Tune WINDOW_SIZE

#define GET_DURATION(START, END) (std::chrono::duration_cast<std::chrono::nanoseconds>((END) - (START)).count() * 1e-9)
#define GET_CURR_TIME() (std::chrono::high_resolution_clock::now())

#define PRUNING_THRESHOLD 10
#define WINDOW_SIZE 1000

LSH::LSH(Dataset &input, const enum metrics &_metric, uint32_t num_ht, uint32_t num_hfs, double _radius):
        maps(new vector< hashtable *>()),
        radius(_radius),
        metric_id(_metric),
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

// Set the appropriate metrics and preprocess the input data
void LSH::set_metrics_and_preprocess() {

    double grid_interval = estimate_grid_interval(raw_inputs);

    if(metric_id == EUCLIDEAN)
        metric = Metrics::Euclidean::distance;
    else {
        for (int i = 0; i < maps->size(); ++i)
            grids->push_back(new Grid(grid_interval));

        // TODO: WATCH OUT THE CAST! IT IS A POSSIBLE FUTURE SEG
        // (if you want your cont. frechet to receive flattened curves as args).
        metric = metric_id == CONTINUOUS_FRECHET ? (distance_f)Metrics::Continuous_Frechet::distance :
                 Metrics::Discrete_Frechet::distance;
    }

    curves_preprocess(*raw_inputs, "input");
}

void LSH::curves_preprocess(std::vector<Curve *> &curves, const std::string& type) {
    for(auto & curve: curves)
        curve_preprocess(*curve, type);
}

template<typename _curve_T>
void LSH::curve_preprocess(_curve_T &c, const string &type) {

    // If the curve already exists, do not re-insert it.
    if(label_to_curve->find(c.get_id()) != label_to_curve->end())
        return;

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
    auto _curve = curve;
    (*grids)[index]->fit(_curve);
    (*grids)[index]->remove_consecutive_duplicates(_curve);
    auto flattened_curve = _curve.flatten();
    flattened_curve->apply_padding(padding_len);
    return flattened_curve;
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

    label_to_curve->clear();
    delete label_to_curve;

    label_to_index_in_fl_queries->clear();
    delete label_to_index_in_fl_queries;
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


    /*double avg_lsh_time_taken = 0;
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
        bruteforce_nn(*(*raw_queries)[i], raw_inputs, Metrics::Discrete_Frechet::distance, &top_brutef);
        end = GET_CURR_TIME();
        auto brutef_time_taken = GET_DURATION(start, end);

        avg_lsh_time_taken += (lsh_time_taken / raw_queries->size());
        abg_brutef_time_taken += (brutef_time_taken / raw_queries->size());

        double candidate_maf = std::get<0>(top_lsh) / std::get<0>(top_brutef);
        if(std::get<0>(top_lsh) < std::numeric_limits<double>::max() - 10e5 && maf < candidate_maf)
            maf = candidate_maf;

        write_data_to_out_file(label, top_lsh, top_brutef, out_path);

        i++;
    }
    write_data_to_out_file(avg_lsh_time_taken, abg_brutef_time_taken, maf, out_path);*/
}

template<typename _curve_T>
void LSH::range_search(_curve_T *query, list<tuple<Curve *, double>> &results) {

    curve_preprocess(*query, "query");
    auto label = query->get_id();
    // Get the query's flattened family from L_flattened_queries
    auto query_family = get_flattened_family(label);
    _range_search(*query_family, results);
}

// Run Range-search
void LSH::_range_search(flattened_curves &query_family, list<tuple<Curve *, double>> &results) {

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

            /* Check for distances, iff items key is equal to query's amplified hf value */
            if(id == identity) {
                // Compute the distance between the query and the current point
                double dist;
                if (dist_cache.find(label) == dist_cache.end()) {

                    Curve *raw_current_curve = label_to_curve->find(label)->second;

                    if(metric_id != EUCLIDEAN) {
                        Curve *raw_query_curve = label_to_curve->find((*query_family[i]).get_id())->second;
                        dist = metric(*raw_query_curve, *raw_current_curve);

                    } else
                        dist = Metrics::Euclidean::distance(*query_family[i], *p);

                    dist_cache.insert({label, dist});

                    if (dist < this->radius)
                        results.emplace_back(raw_current_curve, dist);
                }
            }
        }
    }
}
