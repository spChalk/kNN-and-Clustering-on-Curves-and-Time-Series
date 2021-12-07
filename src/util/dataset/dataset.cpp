
#include <iostream>
#include "dataset.hpp"
#include "../files/file_reader.hpp"

using std::get;
using std::cout;

Dataset::Dataset(std::string & data_path)
: data(file_reader(data_path)) {}

Dataset::~Dataset() {
    for(auto & i : *data)
        delete i;
    delete data;
}

std::vector<Point *> *Dataset::getData() const {
    return data;
}

uint32_t Dataset::get_data_dimensions() {
  auto first_point = (*data)[0];
  return first_point->get_dimensions();
}

uint32_t Dataset::size() {
    return data->size();
}
