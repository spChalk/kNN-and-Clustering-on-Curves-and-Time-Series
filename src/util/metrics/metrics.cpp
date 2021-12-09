
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

// TODO: well, gotta test this...
double Metrics::frechet_distance(Curve &c1, Curve &c2)
{
    auto &a = *c1.get_points();
    auto &b = *c2.get_points();
    uint32_t m1 = a.size();
    uint32_t m2 = b.size();
    double *opt = new double[m1 * m2];

    // Calculate opt for 1st row (index: 0)
    opt[0] = Metrics::euclidean(*a[0], *b[0]);
    auto p1 = a[0];
    for (uint32_t col=1; col != m2; ++col) {
        opt[col] = std::max( Metrics::euclidean(*p1, *b[col]), opt[col-1] );
    }

    // Calculate opt for 1st col (index: row * m2)
    auto q1 = b[0];
    for (uint32_t row=1; row != m1; ++row) {
        uint32_t index = row * m2;
        opt[index] = std::max( Metrics::euclidean(*a[row], *q1), opt[index - m2] );
    }

    // Calculate for rows [2, m1]
    for (uint32_t row=1; row != m1; ++row)
    {
        for (uint32_t col=1; col != m2; ++col)
        {
            auto tmp_min = std::min( opt[(row-1)*m2 + col], std::min( opt[(row-1)*m2 + col - 1], opt[row*m2 + col - 1]) );  // wheelchair
            opt[row * m2 + col] = std::max( Metrics::euclidean(*a[row], *b[col]), tmp_min );
        }
    }

    double ret_val = opt[m1 * m2 - 1];
    delete[] opt;
    return ret_val;
}
