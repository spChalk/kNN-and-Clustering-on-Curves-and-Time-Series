
#ifndef PROJECT_1_METRICS_HPP
#define PROJECT_1_METRICS_HPP

#include "../../curve/point.hpp"
#include "../../curve/curve.hpp"

enum metrics {
    EUCLIDEAN,
    DISCRETE_FRECHET,
    CONTINUOUS_FRECHET
};

class Metrics {
public:

    class Euclidean {
    public:
        static double distance(Point &a, Point &b);
        static double distance(FlattenedCurve &a, FlattenedCurve &b);
        static double distance(Curve &a, Curve &b);
    };

    class Continuous_Frechet {
    public:
        static double distance(Curve &c1, Curve &c2);
        static double distance(FlattenedCurve &c1, FlattenedCurve &c2);
    };

    class Discrete_Frechet {
    public:
        static double distance(Curve &c1, Curve &c2);
    };
};

#endif //PROJECT_1_METRICS_HPP
