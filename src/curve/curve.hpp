#ifndef PROJECT_1_CURVE_HPP
#define PROJECT_1_CURVE_HPP

#include <tuple>
#include <iostream>
#include "../util/metrics/metrics.hpp"

class Curve {

public:
    // TODO DO A DESTRUCTOR
    Curve(std::string id, std::vector<Point *> *points);

    std::string &get_id() {
        return this->id;
    }

    void filter(double pruning_threshold=0.02);
    void erase_time_axis();

    // TODO REMOVE WHEN READY
    void print();

    void fit_to_grid(uint32_t grid_interval);
    std::vector<Point *> *get_points();
    ~Curve();

private:
    std::string id;
    // Data will consist of a vector consisting of points
    std::vector<Point *> *points;

    bool there_is_available_pruning_on_index(uint32_t index, double pruning_threshold);
    void prune_point_on_index(uint32_t index);
    void remove_first_value_of(Point *_tuple);
    uint32_t get_data_dimensions();
};

#endif //PROJECT_1_CURVE_HPP
