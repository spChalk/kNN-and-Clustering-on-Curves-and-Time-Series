
#ifndef PROJECT_1_LSH_HPP
#define PROJECT_1_LSH_HPP

#include "hash_functions/amplified_hf.hpp"
#include "../lsh/hashtable/hashtable.hpp"
#include "../curve/grid.hpp"
#include "../util/metrics/metrics.hpp"
#include <unordered_map>
#include <map>
#include <list>

using std::unordered_multimap;
using std::tuple;
using std::string;
using std::multimap;
using std::vector;
using std::list;

typedef double(*distance_f)(Curve &, Curve &);

// Locality Sensitive Hashing
class LSH {

private:
    // Hash tables
    vector< hashtable *> *maps;

    // Distance function
    distance_f metric;
    const enum metrics metric_id;

    // Radius (for range search)
    double radius;

    std::vector<Grid *> *grids;

    // We keep initial inputs here
    std::vector<Curve *> *raw_inputs;
    std::vector<Curve *> *raw_queries;

    // We store flattened vectors (for L grids each) here
    std::vector< std::vector<FlattenedCurve *> *> *L_flattened_inputs;
    std::vector< std::vector<FlattenedCurve *> *> *L_flattened_queries;

    std::unordered_map<std::string, Curve *> *label_to_curve;

    void curves_preprocess(std::vector<Curve *> &curves, double pr_t, uint32_t max_c_len, const std::string& type);
    void curves_preprocess(std::vector<Curve *> &curves, uint32_t max_c_len, const std::string& type);
    void curves_preprocess(std::vector<Curve *> &curves, const std::string& type, uint32_t num_ht);

    // Run K Nearest Neighbours
    void nn(vector<FlattenedCurve *> &query_family, std::tuple<double, string>& result);

    uint32_t set_metrics_and_preprocess(uint32_t num_ht);
    void delete_flattened_inputs();

public:

    LSH(Dataset &input, Dataset &queries,
        const metrics& metric, uint32_t num_ht = 5, uint32_t num_hfs = 4,
        double radius = 10000);

    // Load data into the structure
    void load();

    void nearest_neighbor(const std::string &out_path);

    // Run Range-search
    void range_search(FlattenedCurve *query, list<tuple<FlattenedCurve *, double>> &);

    virtual ~LSH();

};

#endif //PROJECT_1_LSH_HPP
