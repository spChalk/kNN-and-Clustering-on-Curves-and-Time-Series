
#include <iostream>
#include "dataset.hpp"
#include "../files/file_reader.hpp"

using std::get;
using std::cout;

Dataset::Dataset(std::string & data_path)
: curves(file_reader(data_path)) {}

Dataset::~Dataset() {
    for(auto & i : *curves)
        delete i;
    delete curves;
}

std::vector<Curve *> *Dataset::getData() const {
    return curves;
}

uint32_t Dataset::size() {
    return curves->size();
}

void Dataset::print() {
    for(auto & curve: *curves)
        curve->print();
}

std::vector<FlattenedCurve *> *Dataset::flatten_data() {
    auto fl_data = new std::vector<FlattenedCurve *>();
    for(auto &curve: *getData())
        fl_data->push_back(curve->flatten());
    return fl_data;
}
