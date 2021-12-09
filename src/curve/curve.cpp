#include "../util/metrics/metrics.hpp"
#include "curve.hpp"
#include "grid.hpp"

using std::get;

Curve::Curve(std::string _id, std::vector<Point *> *_points):
id(_id), points(_points) {}

std::string &Curve::get_id() {
    return this->id;
}

void Curve::filter(double pruning_threshold) {

    for(uint32_t i = 1; i < points->size(); i+=2)
        while(there_is_available_pruning_on_index(i, pruning_threshold))
            prune_point_on_index(i);
}

bool Curve::there_is_available_pruning_on_index(uint32_t index, double pruning_threshold) {
    return index+1 < points->size() &&
           Metrics::euclidean(*(*points)[index-1], *(*points)[index]) <= pruning_threshold &&
           Metrics::euclidean(*(*points)[index], *(*points)[index+1]) <= pruning_threshold;
}

void Curve::prune_point_on_index(uint32_t index) {
    delete (*points)[index];
    points->erase(points->begin() + index);
}

void Curve::erase_time_axis() {
    for(auto _tuple: *points)
        remove_first_value_of(_tuple);
}

void Curve::remove_first_value_of(Point *_tuple) {
    auto coordinates = _tuple->get_coordinates();
    coordinates->erase(coordinates->begin());
}

void Curve::print() {
    std::cout << "Curve: " << id << " | ";
    for(auto point: *points)
        point->print();
}

void Curve::fit_to_grid(uint32_t grid_interval) {
    Grid(grid_interval, get_data_dimensions()).fit(*this);
}
/*
vector<Curve *> *Curve::generate_grid_family(uint32_t grid_interval, uint32_t grid_number) {

    auto _curve_on_grids = new vector<Curve *>();
    // First create N grid copies and then populate the
    // "curve_on_grids" class structure
    // in order to avoid eternal loop
    for (int i = 0; i < grid_number; i++) {
        Curve *newCurve = new Curve(*this);
        Grid(grid_interval, get_data_dimensions()).fit(*newCurve);
        _curve_on_grids->push_back(newCurve);
    }
    return _curve_on_grids;
}*/

Curve::~Curve() {
    for(auto & point: *points)
        delete point;
    delete points;
}

uint32_t Curve::get_data_dimensions() {
    return points ? (*points)[0]->get_dimensions() : 0;
}

std::vector<Point *> *Curve::get_points() {
    return points;
}

bool Curve::is_between_min_and_max(uint32_t index) {
    if(index+1 >= points->size())
        return false;

    auto point_left = get_coordinates_of_point(index - 1)[0];
    auto point_middle = get_coordinates_of_point(index)[0];
    auto point_right = get_coordinates_of_point(index + 1)[0];

    if(point_middle >= point_left && point_middle <= point_right ||
        point_middle <= point_left && point_middle >= point_right)
        return true;
    return false;
}

void Curve::min_max_filter() {
    for(uint32_t i = 1; i < points->size(); i+=2)
        while(is_between_min_and_max(i))
            prune_point_on_index(i);
}

std::vector<double> *Curve::get_coordinates_of_point(uint32_t index) {
    assert (index < points->size());
    return (*points)[index]->get_coordinates();
}

void Curve::apply_padding(uint32_t limit) {
    while(points->size() < limit) {
        auto zero = vector<double>{0};
        points->push_back(new Point(zero));
    }
}
