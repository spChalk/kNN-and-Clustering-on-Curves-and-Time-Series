
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
    LSH(Dataset &input, const metrics& metric, uint32_t num_ht = 5, uint32_t num_hfs = 4,
        double radius = 10000);
    LSH(std::vector<FlattenedCurve *> *data, const enum metrics &_metric, uint32_t num_ht = 5, uint32_t num_hfs = 4,
         uint32_t _k_nearest = 1, double _radius = 10000);

    // Load data into the structure
    void load();

    void nearest_neighbor(Curve *query, std::tuple<double, string> &result);

    void range_search(FlattenedCurve *query, std::list<std::tuple<FlattenedCurve *, double>> & top_n);
    void range_search(Curve *query, std::list<std::tuple<Curve *, double>> & top_n);

    virtual ~LSH();

private:
    // Hash tables
    vector< hashtable *> *maps;

    // Distance function
    distance_f metric;
    const enum metrics metric_id;

    // Radius (for range search)
    double radius;

    std::vector<Grid *> *grids = nullptr;

    // We keep initial inputs here
    std::vector<Curve *> *raw_inputs;
    std::vector<Curve *> *raw_queries;
    uint32_t padding_len;

    // We store flattened vectors (for L grids each) here
    vector_of_flattened_curves *L_flattened_inputs = nullptr;
    vector_of_flattened_curves *L_flattened_queries = nullptr;

    std::unordered_map<std::string, Curve *> *label_to_curve = nullptr;
    std::unordered_map<std::string, uint32_t> *label_to_index_in_fl_queries = nullptr;

    flattened_curves *get_flattened_family(std::string &label);

    void curves_preprocess(std::vector<Curve *> &curves, const std::string& type);

    template<typename _curve_T>
    void curve_preprocess(_curve_T &c, const string &type);

    FlattenedCurve *euclidean_preprocess(Curve &curve);
    FlattenedCurve *euclidean_preprocess(FlattenedCurve &curve);
    FlattenedCurve *cont_frechet_preprocess(Curve &curve, uint32_t index);
    FlattenedCurve *discr_frechet_preprocess(Curve &curve, uint32_t index);

    void nn(flattened_curves &query_family, std::tuple<double, string>& result);
    void _range_search(flattened_curves &query_family, std::list<tuple<Curve *, double>> &results);
    void _range_search(FlattenedCurve *query, std::list<std::tuple<FlattenedCurve *, double>> & results);

    template<typename _curve_T>
    void __range_search(_curve_T *query, std::list<tuple<_curve_T *, double>> &results);

    void set_metrics_and_preprocess();

    void delete_flattened_inputs();
};

#endif //PROJECT_1_LSH_HPP
