
#ifndef PROJECT_1_CLUSTER_HPP
#define PROJECT_1_CLUSTER_HPP

#include <chrono>
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <cfloat>
#include <iostream>
#include <iosfwd>
#include <fstream>
#include <iomanip>
#include "../util/dataset/dataset.hpp"
#include "../lsh/lsh.hpp"
#include "../hypercube/hypercube.hpp"
#include "../util/utilities.hpp"
#include "../util/files/file_reader.hpp"

using std::set;
using std::vector;
using std::unordered_map;

class cluster {

private:
    // Pointer to assignment function
    typedef void (cluster::*__CLUSTER_MODULE_assign_func)();
    
    Dataset &dataset;
    /*
     * 0: Classic
     * 1: LSH
     * 2: Hypercube
     */
    uint8_t method;
    uint32_t num_of_clusters;
    double time_taken;

    // Distance function
    distance_f dist_func;
    set< Point *> *centroids;
    set< Point *> *prev_centroids = nullptr;

    LSH *lsh_ds = nullptr;
    hypercube *hc_ds = nullptr;

    double get_min_dist_between_2_centroids();
    double max_distance(unordered_map<Point *, double> &distances);
    void normalize_distances(unordered_map<Point *, double> &distances);
    double compute_sum_of_squared_distances(unordered_map<Point *, double> &distances);
    void assign_point_to_centroid(Point *pt, Point *ctr, std::unordered_map<Point *, std::vector<Point *> *> *final_assign);

    void initialize_centroids();

    std::unordered_map<Point *, std::vector<Point *> *> *final_assign = nullptr;

    void internal_run(__CLUSTER_MODULE_assign_func);

    Point *get_2nd_closest_centroid(Point *, Point *);

    template <typename C>
    void reverse_assignment(C cont);

    void assign_exact_lloyds();
    void assign_LSH();
    void assign_Hypercube();

    double evaluate(std::list<double> &result_per_cluster);

    double update();
    double min_distance_from_centroids(Point *curr_point,
                                       unordered_multimap<string, tuple<string, double>>&);
    void delete_centroids(set< Point *> *);
    void delete_assignment();

public:
    cluster(string &config_path, Dataset &dataset, distance_f, const string &mode="Classic");
    void perform_clustering();
    virtual ~cluster();
    void write_results_to_file(const std::string & out_path, bool verbose=false, bool evalution_on=true);
};

#endif //PROJECT_1_CLUSTER_HPP
