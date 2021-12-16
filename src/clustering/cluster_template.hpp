
#ifndef PROJECT_1_CLUSTER_TMPLT_HPP
#define PROJECT_1_CLUSTER_TMPLT_HPP

#include <chrono>
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <cfloat>
#include <unordered_map>
#include <iostream>
#include <iosfwd>
#include <fstream>
#include <iomanip>
#include "../util/dataset/dataset.hpp"
// #include "../lsh/lsh.hpp"
// #include "../hypercube/hypercube.hpp"
#include "../util/utilities.hpp"
#include "../util/files/file_reader.hpp"
class LSH;
class hypercube; // TODO : RM diz, testing

using std::set;
using std::vector;
using std::unordered_map;

template <typename ItemType, typename DistFunc>
class internal_cluster {

private:
    // Pointer to assignment function
    typedef void (internal_cluster::*__CLUSTER_TMPL_MODULE_assign_func)();
    typedef void (internal_cluster::*__CLUSTER_TMPL_MODULE_update_func)();
    
    Dataset &dataset;

    uint8_t assign_method;  /* 0: Classic, 1: LSH, 2: Hypercube */
    uint8_t update_method;  /* 0: Vector,  1: Curve             */

    uint32_t num_of_clusters;
    double time_taken;

    // Distance function
    DistFunc dist_func;
    set< ItemType *> *centroids;
    set< ItemType *> *prev_centroids = nullptr;

    LSH *lsh_ds = nullptr;
    hypercube *hc_ds = nullptr;

    double get_min_dist_between_2_centroids();
    double max_distance(unordered_map<ItemType *, double> &distances);
    void normalize_distances(unordered_map<ItemType *, double> &distances);
    double compute_sum_of_squared_distances(unordered_map<ItemType *, double> &distances);
    void assign_point_to_centroid(ItemType *pt, ItemType *ctr, std::unordered_map<ItemType *, std::vector<ItemType *> *> *final_assign);

    void initialize_centroids();

    std::unordered_map<ItemType *, std::vector<ItemType *> *> *final_assign = nullptr;

    void internal_run(__CLUSTER_TMPL_MODULE_assign_func, __CLUSTER_TMPL_MODULE_update_func);

    ItemType *get_2nd_closest_centroid(ItemType *, ItemType *);

    template <typename C>
    void reverse_assignment(C cont);

    void assign_exact_lloyds();
    void assign_LSH();
    void assign_Hypercube();

    double evaluate(std::list<double> &result_per_cluster);

    double update_vector();
    double update_curve();
    double min_distance_from_centroids(ItemType *curr_point,
                                       std::unordered_multimap<string, std::tuple<string, double>>&);
    void delete_centroids(set< ItemType *> *);
    void delete_assignment();

public:
    internal_cluster(string &config_path, Dataset &dataset, DistFunc dist_func, const string &mode="Classic");
    ~internal_cluster();
    void perform_clustering();
    void write_results_to_file(const std::string & out_path, bool verbose=false, bool evalution_on=true);
};

#endif //PROJECT_1_CLUSTER_TMPLT_HPP
