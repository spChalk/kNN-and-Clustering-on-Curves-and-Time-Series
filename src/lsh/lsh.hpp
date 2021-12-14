
#ifndef PROJECT_1_LSH_HPP
#define PROJECT_1_LSH_HPP

#include "hash_functions/amplified_hf.hpp"
#include "../lsh/hashtable/hashtable.hpp"
#include "../curve/grid.hpp"
#include <unordered_map>
#include <map>
#include <list>

using std::unordered_multimap;
using std::tuple;
using std::string;
using std::multimap;
using std::vector;
using std::list;

typedef double(*distance_f)(FlattenedCurve&, FlattenedCurve&);

// Locality Sensitive Hashing
class LSH {

private:
    // Hash tables
    vector< hashtable *> *maps;
    // Distance function
    distance_f distance;
    // Radius (for range search)
    double radius;

    Grid *grid;

    std::vector<FlattenedCurve *> *flattened_input_data;
    std::vector<FlattenedCurve *> *flattened_query_data;

    void curves_preprocess(std::vector<Curve *> &curves, double pruning_threshold, uint32_t max_curve_length);

    // Run K Nearest Neighbours
    void nn(FlattenedCurve &query, std::tuple<double, string>& result);

public:

    LSH(Dataset &input, Dataset &queries,
        distance_f, uint32_t num_ht = 5, uint32_t num_hfs = 4,
        double radius = 10000);

    // Load data into the structure
    void load(vector<FlattenedCurve *> *data);

    void nearest_neighbor(const std::string &out_path);

    // Run Range-search
    void range_search(FlattenedCurve *query, list<tuple<FlattenedCurve *, double>> &);

    virtual ~LSH();
};

#endif //PROJECT_1_LSH_HPP
