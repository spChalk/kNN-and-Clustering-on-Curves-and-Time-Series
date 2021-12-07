#ifndef PROJECT_1_CURVE_HPP
#define PROJECT_1_CURVE_HPP

#include <tuple>
#include <iostream>
#include "../util/dataset/dataset.hpp"
#include "../util/metrics/metrics.hpp"

class Curve: public Dataset {

public:
    Curve(std::string & data_path);

    void filter(double pruning_threshold=0.02);
    void erase_time_axis();
    // TODO REMOVE WHEN READY
    void print();

    void fit_to_grids(uint32_t grid_interval, uint32_t grid_number);

    ~Curve();
private:
    // Grid family
    std::vector<Curve *> *curve_on_grids;

    bool there_is_available_pruning_on_index(uint32_t index, double pruning_threshold);
    void prune_point_on_index(uint32_t index);
    void remove_first_value_of(Point *_tuple);
    std::vector<Curve *> *generate_grid_family(uint32_t grid_interval, uint32_t grid_number);

};

#endif //PROJECT_1_CURVE_HPP
