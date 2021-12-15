
#ifndef PROJECT_1_LSH_HPP
#define PROJECT_1_LSH_HPP

#include "hash_functions/amplified_hf.hpp"
#include "../lsh/hashtable/hashtable.hpp"
#include "../curve/grid.hpp"
#include "../util/metrics/metrics.hpp"
#include <unordered_map>
#include <map>
#include <list>

typedef double(*distance_f)(Curve &, Curve &);
typedef std::vector<FlattenedCurve *> flattened_curves;
typedef std::vector< flattened_curves *> vector_of_flattened_curves;

// Locality Sensitive Hashing
class LSH {

public:
    LSH(Dataset &input, Dataset &queries, const metrics& metric, uint32_t num_ht = 5, uint32_t num_hfs = 4,
        double radius = 10000);

    // Load data into the structure
    void load();

    void nearest_neighbor(const std::string &out_path);

    // Run Range-search
    void range_search(flattened_curves &query_family, std::list<tuple<Curve *, double>> &results);

    virtual ~LSH();

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
    vector_of_flattened_curves *L_flattened_inputs;
    vector_of_flattened_curves *L_flattened_queries;

    std::unordered_map<std::string, Curve *> *label_to_curve;

    void curves_preprocess(std::vector<Curve *> &curves, double pr_t, uint32_t max_c_len, const std::string& type);
    void curves_preprocess(std::vector<Curve *> &curves, uint32_t max_c_len, const std::string& type);
    void curves_preprocess(std::vector<Curve *> &curves, const std::string& type, uint32_t num_ht);

    // Run K Nearest Neighbours
    void nn(flattened_curves &query_family, std::tuple<double, string>& result);

    uint32_t set_metrics_and_preprocess(uint32_t num_ht);

    void delete_flattened_inputs();
};

#endif //PROJECT_1_LSH_HPP
