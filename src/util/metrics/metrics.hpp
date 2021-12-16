
#ifndef PROJECT_1_METRICS_HPP
#define PROJECT_1_METRICS_HPP

#include "../../curve/point.hpp"
#include "../../curve/curve.hpp"

#include <list>
#include <tuple>

enum metrics {
    EUCLIDEAN,
    DISCRETE_FRECHET,
    CONTINUOUS_FRECHET
};

namespace Metrics {

    namespace Euclidean {
        double distance(Point &a, Point &b);
        double distance(FlattenedCurve &a, FlattenedCurve &b);
        double distance(Curve &a, Curve &b);
    };

    namespace Continuous_Frechet {
        double distance(Curve &c1, Curve &c2);
        double distance(FlattenedCurve &c1, FlattenedCurve &c2);
    };

    namespace Discrete_Frechet {
        double distance(Curve &c1, Curve &c2);
        void optimal_traversal(Curve &c1, Curve &c2, std::list<std::tuple<uint32_t, uint32_t>> &lp);
        void clean();
    };
};

#endif //PROJECT_1_METRICS_HPP
