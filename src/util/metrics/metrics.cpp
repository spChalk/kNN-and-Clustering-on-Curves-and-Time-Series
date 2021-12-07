
#include "metrics.hpp"

#include <iostream>
#include <iosfwd>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <cmath>

double Metrics::euclidean(Point &a, Point &b) {

    if(a.get_dimensions() != b.get_dimensions()) {
        std::ostringstream msg;
        msg << "Exception in 'distance'. Vectors do not have the same size." << std::endl;
        throw std::runtime_error(msg.str());
    }
    auto a_coord = *a.get_coordinates();
    auto b_coord = *b.get_coordinates();
    double result = 0;
    for(uint32_t i = 0; i < a.get_dimensions(); i++) {
        double factor = a_coord[i] - b_coord[i];
        result += (factor * factor);
    }
    return sqrt(result);
}