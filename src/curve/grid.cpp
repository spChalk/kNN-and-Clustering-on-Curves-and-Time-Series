#include <vector>
#include "grid.hpp"

using std::vector;
using std::get;

Grid::Grid(double _grid_interval, uint32_t dimensions):
grid_interval(_grid_interval), noise(Distributions::uniform<double>(0, (long)_grid_interval)) {
   //initialize_noise(dimensions);
}

void Grid::fit(Curve &curve) {
    map_to_grid(curve);
    //remove_consecutive_duplicates(curve);
    //add_noise(curve);
}
/*
void Grid::initialize_noise(uint32_t dimensions) {
    for(uint32_t i = 0; i < dimensions; i++)
        noise->push_back(Distributions::uniform<double>(0, this->grid_interval));
}

void Grid::add_noise(Curve &curve) {
    for(auto & point: *curve.getData())
        add_vectors<double>(point->get_coordinates(), noise, point->get_coordinates());
}

void Grid::remove_consecutive_duplicates(Curve &curve) {
    auto curve_data = curve.getData();
    for(uint32_t i=0; i < curve.size() - 1; i++)
        while (are_equal_consecutive_vectors(curve_data, i))
            curve_data->erase(curve_data->begin() + i + 1);
}

bool Grid::are_equal_consecutive_vectors(const std::vector<Point *> *curve_data,
                                        uint32_t index) const {
    return vectors_are_equal<double>((*curve_data)[index]->get_coordinates(), (*curve_data)[index + 1]->get_coordinates());
}
*/
void Grid::snap(vector<double> *_vector, double interval) {
    // TODO: Maybe overload operators in order to simplify this circus below
    subtract_scalar_from_vector<double>(_vector, noise);
    divide_vector_by_scalar<double>(_vector, interval);
    add_scalar_to_vector<double>(_vector, 0.5);
    apply_floor_to_vector(_vector);
    multiply_vector_by_scalar<double>(_vector, interval);
    add_scalar_to_vector<double>(_vector, noise);
}

void Grid::map_to_grid(Curve &curve) {
    for(auto &_point: *curve.get_points())
        snap(_point->get_coordinates(), grid_interval);
}
