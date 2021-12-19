#include <iostream>

#include "acutest.h"
#include "../src/util/metrics/metrics.hpp"


void test_dfd_and_traversal(void);
Curve *create_curve(int num_points);
double check_traversal(Curve &c1, Curve &c2, std::list<std::tuple<uint32_t, uint32_t>> &lp);
double get_dfd_bruteforce(Curve &c1, uint32_t index1, Curve &c2, uint32_t index2);
void test_pair_curves(int dim1, int dim2);


TEST_LIST = {
    { "test_dfd", test_dfd_and_traversal },
    { NULL, NULL }
};

void test_dfd_and_traversal(void)
{
    int num_tests = 50;
    int max_dim = 10;
    int min_dim = 2;
    for (int i = 0; i < num_tests; ++i)
    {
        int dim1 = rand() % max_dim + min_dim;
        int dim2 = rand() % max_dim + min_dim;
        test_pair_curves(dim1, dim2);
    }
}

void test_pair_curves(int dim1, int dim2)
{
    Curve *c1 = create_curve(dim1);
    Curve *c2 = create_curve(dim2);

    double dist_bf = get_dfd_bruteforce(*c1, 0, *c2, 0);
    double dist_dfd = Metrics::Discrete_Frechet::distance(*c1, *c2);

    auto lp = std::list<std::tuple<uint32_t, uint32_t>>();
    Metrics::Discrete_Frechet::optimal_traversal(*c1, *c2, lp);
    double dfd_trav = check_traversal(*c1, *c2, lp);

    TEST_ASSERT(dist_bf == dist_dfd);
    TEST_ASSERT(dist_bf == dfd_trav);

    Metrics::Discrete_Frechet::clean();
    delete c1;
    delete c2;
}

Curve *create_curve(int num_points) {
    static int id_counter = 0;
    auto *vec_with_points = new std::vector<Point *>();
    for (int i = 0; i < num_points; i++) {
        double rand_num = rand() % 10 + 1;
        double x_item = i+1;
        auto vec = std::vector<double>{x_item, rand_num};
        Point *p = new Point(vec);
        vec_with_points->emplace_back(p);
    }
    auto name = std::to_string(id_counter++);
    Curve *c = new Curve(name, vec_with_points);
    return c;
}


double get_dfd_bruteforce(Curve &c1, uint32_t index1, Curve &c2, uint32_t index2) {
    // 3 options:
    // move C1, move C2, move both by 1
    auto p1 = c1.get_points()->at(index1);
    auto q1 = c2.get_points()->at(index2);

    double curr_dist = Metrics::Euclidean::distance(*p1, *q1);

    if (index1 == c1.get_points()->size()-1 && index2 == c2.get_points()->size()-1) {
        return curr_dist;
    }
    if (index1 == c1.get_points()->size()-1) {
        return std::max(curr_dist, get_dfd_bruteforce(c1, index1, c2, index2+1));
    }
    if (index2 == c2.get_points()->size()-1) {
        return std::max(curr_dist, get_dfd_bruteforce(c1, index1+1, c2, index2));
    }
   
    double moveC1 = get_dfd_bruteforce(c1, index1+1, c2, index2);
    double moveC2 = get_dfd_bruteforce(c1, index1, c2, index2+1);
    double moveBT = get_dfd_bruteforce(c1, index1+1, c2, index2+1);

    double min1 = std::min(moveC1, moveC2);
    double final_min = std::min(min1, moveBT);
    return std::max(final_min, curr_dist);
}

double check_traversal(Curve &c1, Curve &c2, std::list<std::tuple<uint32_t, uint32_t>> &lp)
{
    double max_pair_dist = -1.0;
    uint32_t index1 = (uint32_t)(-1);
    uint32_t index2 = (uint32_t)(-1);

    for (auto it = lp.begin(); it != lp.end(); ++it)
    {
        auto &t = *it;
        uint32_t x = std::get<0>(t);
        uint32_t y = std::get<1>(t);
        
        TEST_ASSERT(x == index1 || x == index1 + 1);
        TEST_ASSERT(y == index2 || y == index2 + 1);
        TEST_ASSERT(!(x == index1 && y == index2));

        auto *p = c1.get_points()->at(x);
        auto *q = c2.get_points()->at(y);
        double dist = Metrics::Euclidean::distance(*p, *q);
        if (dist > max_pair_dist)
            max_pair_dist = dist;
        
        index1 = x;
        index2 = y;
    }

    TEST_ASSERT(index1 == c1.get_points()->size()-1);
    TEST_ASSERT(index2 == c2.get_points()->size()-1);
    return max_pair_dist;   
}