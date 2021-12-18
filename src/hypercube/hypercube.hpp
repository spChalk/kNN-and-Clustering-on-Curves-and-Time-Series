
#ifndef PROJECT_1_HYPERCUBE_HPP
#define PROJECT_1_HYPERCUBE_HPP

#include <list>
#include <vector>
#include <tuple>
#include <queue>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "../util/dataset/dataset.hpp"
#include "../lsh/hash_functions/hash_function.hpp"

// The Hypercube hash table is a vector, and each bucket is also represented as a vector
// (seperate-chaining HT implementation)
typedef std::vector<std::vector<FlattenedCurve *>> hash_table_hypercube;

// User-given distance metric
typedef double(*hypercube_dist_func)(FlattenedCurve&, FlattenedCurve&);

class hypercube
{
private:
    // Hypercube (data structure) parameters
    const uint32_t init_dim;  // Initial dimension (d) before projection.
    const uint32_t k;         // k: Dimensions projected (d') (d' << d)
    hypercube_dist_func dist_func;     // Distance metric 

    // Parameters of NN queries
    uint32_t M;       // Max points checked
    uint32_t probes;  // Max vertices checked
    uint32_t N;       // Get N nearest neighbors
    uint32_t R;       // Radius

    hash_table_hypercube *hypercube_ht;
    std::unordered_map<uint64_t, int> *bucket_to_vertex_index;
    std::vector<hash_function *> *hash_family;  // LSH-admitting hash family

    // Performs a bucket projection to a hypercube vertex
    uint32_t convert_bucket_to_vertex(FlattenedCurve *datapoint);
    uint32_t select_next_vertex(bool & visited_all, uint32_t current_vertex, uint32_t max_vertices, std::queue<uint32_t> &next_vertices, std::unordered_set<uint32_t> & vertex_cache);

    template <typename C, typename T>
    void perform_query(FlattenedCurve *query, T metric, hypercube_dist_func dist_func, C & top_n);

public:
    // Structural Hypercube operations
    hypercube(Dataset &dataset, hypercube_dist_func dist_func, uint32_t k=14, uint32_t M=10, uint32_t probes=2, uint32_t N=1, uint32_t R=10000);
    void set_limits(uint32_t N=1, uint32_t R=10000, uint32_t M=10, uint32_t probes=2);
    ~hypercube();

    void insert_point(FlattenedCurve *);
 
    // Query operations
    uint32_t search_for_vertex(FlattenedCurve *p);
    void knn(FlattenedCurve *query, std::multimap<double, std::string> & top_n);
    void range_search(FlattenedCurve *query, std::list<std::tuple<FlattenedCurve *, double>> & top_n);
    void range_search(Curve *query, std::list<std::tuple<Curve *, double>> & top_n){};
};



#endif //PROJECT_1_HYPERCUBE_HPP
