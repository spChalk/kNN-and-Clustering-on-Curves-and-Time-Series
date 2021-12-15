#include <vector>
#include "grid.hpp"

using std::vector;
using std::get;

Grid::Grid(double _grid_interval):
grid_interval(_grid_interval), noise(Distributions::uniform<double>(0, (long)_grid_interval)) {
}

void Grid::fit(Curve &curve) {
    map_to_grid(curve);
}

void Grid::remove_consecutive_duplicates(Curve &curve) {
    auto curve_data = curve.get_points();
    for(uint32_t i=0; i < curve.get_points()->size() - 1; i++)
        while (i+1 < curve.get_points()->size() &&
                    are_equal_consecutive_vectors(curve_data, i)) {
            auto to_delete = (*curve_data)[i+1];
            curve_data->erase(curve_data->begin() + i + 1);
            delete to_delete;
        }
}

bool Grid::are_equal_consecutive_vectors(std::vector<Point *> *curve_data,
                                        uint32_t index) {
    return vectors_are_equal<double>((*curve_data)[index]->get_coordinates(), (*curve_data)[index + 1]->get_coordinates());
}

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
