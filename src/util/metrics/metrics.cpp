
#include "metrics.hpp"
#include "../../../cont_frechet_repo/Fred/include/frechet.hpp"

#include <iostream>
#include <iosfwd>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <memory>

double Metrics::Euclidean::distance(Curve &a, Curve &b) {
    auto ap = a.get_points();
    auto bp = b.get_points();
    uint32_t size = std::min( ap->size(), bp->size() );
    if (size == 0) return 0;
    double total = 0;
    for (uint32_t i = 0; i != size; ++i) {
        total += Metrics::Euclidean::distance(*ap->at(i), *bp->at(i));
    }
    return total / size;
}

double Metrics::Euclidean::distance(FlattenedCurve &a, FlattenedCurve &b) {

    if(a.get_size() != b.get_size()) {
        std::ostringstream msg;
        msg << "Exception in 'distance'. Curves do not have the same size." << std::endl;
        throw std::runtime_error(msg.str());
    }
    auto a_coord = *a.get_coordinates();
    auto b_coord = *b.get_coordinates();
    double result = 0;
    for(uint32_t i = 0; i < a.get_size(); i++) {
        double factor = a_coord[i] - b_coord[i];
        result += (factor * factor);
    }
    return sqrt(result);
}

double Metrics::Euclidean::distance(Point &a, Point &b) {

    if(a.get_dimensions() != b.get_dimensions()) {
        std::ostringstream msg;
        msg << "Exception in 'distance'. Curves do not have the same size." << std::endl;
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

static double *opt_array = nullptr;

// TODO: well, gotta test this...
double Metrics::Discrete_Frechet::distance(Curve &c1, Curve &c2)
{
    auto &a = *c1.get_points();
    auto &b = *c2.get_points();
    uint32_t m1 = a.size();
    uint32_t m2 = b.size();
    

    auto& opt = opt_array;

    if (opt != nullptr)
        delete[] opt;

    opt = new double[m1 * m2];

    std::shared_ptr<double> shared_good(opt, std::default_delete<double[]>());


    // Calculate opt for 1st row (index: 0)
    opt[0] = Metrics::Euclidean::distance(*a[0], *b[0]);
    auto p1 = a[0];
    for (uint32_t col=1; col != m2; ++col) {
        opt[col] = std::max( Metrics::Euclidean::distance(*p1, *b[col]), opt[col-1] );
    }

    // Calculate opt for 1st col (index: row * m2)
    auto q1 = b[0];
    for (uint32_t row=1; row != m1; ++row) {
        uint32_t index = row * m2;
        opt[index] = std::max( Metrics::Euclidean::distance(*a[row], *q1), opt[index - m2] );
    }

    // Calculate for rows [2, m1]
    for (uint32_t row=1; row != m1; ++row)
    {
        for (uint32_t col=1; col != m2; ++col)
        {
            auto tmp_min = std::min( opt[(row-1)*m2 + col], std::min( opt[(row-1)*m2 + col - 1], opt[row*m2 + col - 1]) );  // wheelchair
            opt[row * m2 + col] = std::max( Metrics::Euclidean::distance(*a[row], *b[col]), tmp_min );
        }
    }

    double ret_val = opt[m1 * m2 - 1];
    return ret_val;
}

#define MAX_OF_THREE(X, Y, Z) { \
    X > Y ? (X > Z ? 0 : 2) \
          : (Y > Z ? 1 : 2) \
}

void Metrics::Discrete_Frechet::optimal_traversal(Curve &c1, Curve &c2, std::list<std::tuple<uint32_t, uint32_t>> &lp)
{
    auto& opt = opt_array;

    auto &a = *c1.get_points();
    auto &b = *c2.get_points();
    uint32_t m1 = a.size();
    uint32_t m2 = b.size();

    uint32_t p = m1, q = m2;  // TODO : test this - needs tweaking
    std::tuple<uint32_t, uint32_t> t = std::make_tuple(p, q);
    lp.push_front(t);

    while (p && q)
    {
        // TODO : test this - needs tweaking
        double x = opt[p * m2 + q];
        double y = opt[p * m2 + q-1];
        double z = opt[(p-1)*m2 + q-1];
        uint32_t minIdx = (int) MAX_OF_THREE(x, y, z);

        if (minIdx == 0) {
            t = std::make_tuple(--p, q);
        }
        else if (minIdx == 1) {
            t = std::make_tuple(p, --q);
        }
        else {
            t = std::make_tuple(--p, --q);
        }
        lp.push_front(t);
    }

    delete[] opt;
}

void Metrics::Discrete_Frechet::clean() { delete[] opt_array; }


double Metrics::Continuous_Frechet::distance(Curve &c1, Curve &c2) {
    _Curve *fc1 = c1.to_FredCurve();
    _Curve *fc2 = c2.to_FredCurve();
    auto result = Frechet::Continuous::distance(*fc1, *fc2).value;
    delete fc1; delete fc2;
    return result;
}

double Metrics::Continuous_Frechet::distance(FlattenedCurve &c1, FlattenedCurve &c2) {
    _Curve *fc1 = c1.to_FredCurve();
    _Curve *fc2 = c2.to_FredCurve();
    auto result = Frechet::Continuous::distance(*fc1, *fc2).value;
    delete fc1; delete fc2;
    return result;
}

