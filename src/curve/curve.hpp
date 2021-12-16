#ifndef PROJECT_1_CURVE_HPP
#define PROJECT_1_CURVE_HPP

#include <tuple>
#include <iostream>
#include "./point.hpp"
#include "../../cont_frechet_repo/Fred/include/curve.hpp"

class FlattenedCurve;
class Curve {

public:
    // TODO DO A DESTRUCTOR
    Curve(std::string &id, std::vector<Point *> *points);
    Curve(const Curve &curve);

    void filter(double pruning_threshold=0.02);
    void erase_time_axis();
    void min_max_filter();
    void apply_padding(uint32_t limit);

    std::string &get_id();
    std::vector<Point *> *get_points();
    std::vector<double> *get_coordinates_of_point(uint32_t index);

    ~Curve();

    // TODO REMOVE WHEN READY
    void print();
    uint32_t get_data_dimensions();

    _Curve *to_FredCurve();
    FlattenedCurve *flatten();

private:
    std::string id;
    // Data will consist of a vector consisting of points
    std::vector<Point *> *points;

    bool there_is_available_pruning_on_index(uint32_t index, double pruning_threshold);
    void prune_point_on_index(uint32_t index);
    void remove_first_value_of(Point *_tuple);
    bool is_between_min_and_max(uint32_t index);
};

class FlattenedCurve {
private:
    std::string id;
    std::vector<double> *points;

public:
    FlattenedCurve(Curve &normal_curve);
    FlattenedCurve(const FlattenedCurve &curve);
    FlattenedCurve(std::string _id, std::vector<double> &_coord)
            : id(_id), points(new std::vector<double>(_coord)) {}

    ~FlattenedCurve();

    std::string get_id();
    std::vector<double> *get_coordinates();
    uint32_t get_size();
    void apply_padding(uint32_t limit);

    _Curve *to_FredCurve();
};

#endif //PROJECT_1_CURVE_HPP
