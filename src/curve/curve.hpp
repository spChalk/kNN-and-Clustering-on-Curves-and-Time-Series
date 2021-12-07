#ifndef PROJECT_1_CURVE_HPP
#define PROJECT_1_CURVE_HPP

#include <tuple>
#include <iostream>
#include "../util/metrics/metrics.hpp"

class Curve {

public:
    // TODO DO A DESTRUCTOR
    Curve(std::string id, std::vector<Point *> *points);

    void filter(double pruning_threshold=0.02);
    void erase_time_axis();
    void fit_to_grid(uint32_t grid_interval);
    void min_max_filter();
    void apply_padding(uint32_t limit);

    std::string &get_id();
    std::vector<Point *> *get_points();
    std::vector<double> *get_coordinates_of_point(uint32_t index);

    ~Curve();

    // TODO REMOVE WHEN READY
    void print();

private:
    std::string id;
    // Data will consist of a vector consisting of points
    std::vector<Point *> *points;

    bool there_is_available_pruning_on_index(uint32_t index, double pruning_threshold);
    void prune_point_on_index(uint32_t index);
    void remove_first_value_of(Point *_tuple);
    bool is_between_min_and_max(uint32_t index);
    uint32_t get_data_dimensions();
};

#endif //PROJECT_1_CURVE_HPP
