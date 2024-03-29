#include "../util/metrics/metrics.hpp"
#include "curve.hpp"
#include "grid.hpp"

using std::get;

Curve::Curve(std::string &_id, std::vector<Point *> *_points):
id(_id), points(_points), x_axis_erased(false) {}

Curve::Curve(const Curve &curve):
id(curve.id), points(new std::vector<Point *>()), x_axis_erased(curve.x_axis_erased) {

    for(auto &point: *curve.points)
        points->push_back(new Point(*point->get_coordinates()));
}

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
           Metrics::Euclidean::distance(*(*points)[index-1], *(*points)[index]) <= pruning_threshold &&
           Metrics::Euclidean::distance(*(*points)[index], *(*points)[index+1]) <= pruning_threshold;
}

void Curve::prune_point_on_index(uint32_t index) {
    delete (*points)[index];
    points->erase(points->begin() + index);
}

void Curve::erase_time_axis() {
    if (!this->x_axis_erased) {
        this->x_axis_erased = true;
        for(auto _tuple: *points)
            remove_first_value_of(_tuple);
    }
}

void Curve::remove_first_value_of(Point *_tuple) {
    auto coordinates = _tuple->get_coordinates();
    coordinates->erase(coordinates->begin());
}

uint32_t Curve::get_size() {
    return points->size();
}


void Curve::print() {
    std::cout << "Curve: " << id << " |\n";
    for(auto point: *points)
        point->print();
}

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

    if((point_middle >= point_left && point_middle <= point_right) ||
        (point_middle <= point_left && point_middle >= point_right))
        return true;
    return false;
}

void Curve::min_max_filter() {
    for(uint32_t i = 1; i < points->size(); i+=2)
        while(is_between_min_and_max(i))
            prune_point_on_index(i);
}

std::vector<double> *Curve::get_coordinates() {
    // This function exists for template compatibility
    throw std::runtime_error("Curve::get_coordinates() called");
    return this->get_coordinates_of_point(0);
}


std::vector<double> *Curve::get_coordinates_of_point(uint32_t index) {
    assert (index < points->size());
    return (*points)[index]->get_coordinates();
}

void Curve::apply_padding(uint32_t limit) {
    uint32_t point_dim = this->points->at(0)->get_dimensions();
    while(points->size() < limit) {
        double pad_value = 1e9;
        auto pad_vec = vector<double>(point_dim, pad_value);
        points->push_back(new Point(pad_vec));
    }
}

// This function exists for template compatibility
Curve::Curve(std::string &_id, std::vector<double> &first_point)
:id(_id), points(nullptr) { throw std::runtime_error("Curve::get_coordinates() called"); }

_Curve *Curve::to_FredCurve() {

    _Curve *fredCurve = new _Curve(Points(points->size(), get_data_dimensions()));
    for (uint32_t i = 0; i < fredCurve->size(); i++) {
        for (uint32_t j = 0; j < fredCurve->operator[](i).dimensions(); j++) {
            fredCurve->operator[](i).set(j, (*(*points)[i]->get_coordinates())[j]);
        }
    }
    return fredCurve;
}

FlattenedCurve *Curve::flatten() {
    return new FlattenedCurve(*this);
}

FlattenedCurve::FlattenedCurve(Curve &normal_curve)
: id(normal_curve.get_id()), points(new vector<double>()) {
    for(auto &point: *normal_curve.get_points()) {
        for(auto &coord: *point->get_coordinates()) {
            points->push_back(coord);
        }
    }
}

FlattenedCurve::FlattenedCurve(const FlattenedCurve& curve):
id(curve.id), points(new vector<double>()) {
    for(auto &coord: *curve.points) {
        this->points->push_back(coord);
    }
}

FlattenedCurve::~FlattenedCurve() {
    delete points;
}

std::string FlattenedCurve::get_id() {
    return id;
}

std::vector<double> *FlattenedCurve::get_coordinates() {
    return points;
}

uint32_t FlattenedCurve::get_size() {
    return points->size();
}

uint32_t FlattenedCurve::get_data_dimensions() {
    return points->size();
}

_Curve *FlattenedCurve::to_FredCurve() {
    _Curve *fredCurve = new _Curve(Points(points->size(), 1));
    for (uint32_t i = 0; i < fredCurve->size(); i++)
        fredCurve->operator[](i).set(0, (*points)[i]);
    return fredCurve;
}

void FlattenedCurve::apply_padding(uint32_t limit) {
    while(points->size() < limit) {
        points->push_back(0.0);
    }
}
