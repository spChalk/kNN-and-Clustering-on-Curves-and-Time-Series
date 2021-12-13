
#ifndef PROJECT_1_LSH_HPP
#define PROJECT_1_LSH_HPP

#include "hash_functions/amplified_hf.hpp"
#include "../lsh/hashtable/hashtable.hpp"
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
    // Number of nearest neighbours
    uint32_t k_nearest_n;
    // Radius (for range search)
    double radius;

public:

    LSH(vector<FlattenedCurve *> * data, distance_f, uint32_t num_ht = 5,
        uint32_t num_hfs = 4, uint32_t k_nearest = 1, double radius = 10000);

    // Load data into the structure
    void load(vector<FlattenedCurve *> *data);

    // Run K Nearest Neighbours
    void knn(FlattenedCurve *query, multimap<double, string> &);

    // Run Range-search
    void range_search(FlattenedCurve *query, list<tuple<FlattenedCurve *, double>> &);

    virtual ~LSH();
};

#endif //PROJECT_1_LSH_HPP
