
#include "curve.hpp"
#include "grid.hpp"

using std::get;


Curve::Curve(string &data_path):
Dataset(data_path), curve_on_grids(nullptr) {}

void Curve::filter(double pruning_threshold) {

    for(uint32_t i = 1; i < data->size(); i+=2)
        while(there_is_available_pruning_on_index(i, pruning_threshold))
            prune_point_on_index(i);
}

bool Curve::there_is_available_pruning_on_index(uint32_t index, double pruning_threshold) {
    return index+1 < data->size() &&
           Metrics::euclidean(*(*data)[index-1], *(*data)[index]) <= pruning_threshold &&
           Metrics::euclidean(*(*data)[index], *(*data)[index+1]) <= pruning_threshold;
}

void Curve::prune_point_on_index(uint32_t index) {
    delete (*data)[index];
    data->erase(data->begin() + index);
}

void Curve::erase_time_axis() {
    for(auto _tuple: *data)
        remove_first_value_of(_tuple);
}

void Curve::remove_first_value_of(Point *_tuple) {
    auto coordinates = _tuple->get_coordinates();
    coordinates->erase(coordinates->begin());
}

void Curve::print() {

    for(auto t: *data) {
        std::cout << "( ";
        std::cout << t->get_id() << " , ";
        for(auto i: *t->get_coordinates()) {
            std::cout << i << " ";
        }
        std::cout << ")\n";
    }
}

void Curve::fit_to_grids(uint32_t grid_interval, uint32_t grid_number) {

    this->curve_on_grids =generate_grid_family(grid_interval, grid_number);

    // TODO Apply padding to the resulting vectors
}

vector<Curve *> *Curve::generate_grid_family(uint32_t grid_interval, uint32_t grid_number) {

    auto _curve_on_grids = new vector<Curve *>();
    // First create N grid copies and then populate the
    // "curve_on_grids" class structure
    // in order to avoid eternal loop
    for (int i = 0; i < grid_number; i++) {
        //TODO CREATE COPY CONSTRUCTOR
        Curve *newCurve = new Curve(*this);
        Grid(grid_interval, get_data_dimensions()).fit(*newCurve);
        _curve_on_grids->push_back(newCurve);
    }
    return _curve_on_grids;
}

Curve::~Curve() {
    //TODO
}