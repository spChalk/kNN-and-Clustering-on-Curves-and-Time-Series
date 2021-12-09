
#ifndef PROJECT_1_METRICS_HPP
#define PROJECT_1_METRICS_HPP

#include "../../curve/point.hpp"
#include "../../curve/curve.hpp"

class Metrics {

public:
    static double euclidean(Point &a, Point &b);
    static double frechet_distance(Curve &c1, Curve &c2);
    // static cont_frechet
};

#endif //PROJECT_1_METRICS_HPP
