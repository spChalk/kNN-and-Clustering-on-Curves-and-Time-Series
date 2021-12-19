#include <iostream>

#include "acutest.h"
#include "../src/util/metrics/metrics.hpp"
#include "../src/clustering/complete_binary_tree.hpp"

void test_cbt_post_order_sum(void);
int *internal_post_order(CompleteBinaryTree<int> &tree, CBTree_Node node);
int post_order(CompleteBinaryTree<int> &tree);
int inner_test_cbt(int max_value);

void test_dfd(void);
void test_dfd_and_traversal(void);
Curve *create_curve(int num_points);
double check_traversal(Curve &c1, Curve &c2, std::list<std::tuple<uint32_t, uint32_t>> &lp);
double get_dfd_bruteforce(Curve &c1, uint32_t index1, Curve &c2, uint32_t index2);
void test_pair_curves(int dim1, int dim2, bool check_opt_traversal);


TEST_LIST = {
    { "test_cbt", test_cbt_post_order_sum },
    { "test_dfd", test_dfd },
    { "test_opt_traversal", test_dfd_and_traversal },
    { NULL, NULL }
};

/* Discrete Frechet Distance and Optimal Traversal TESTS */

void test_dfd(void)
{
    int num_tests = 50;
    int max_dim = 10;
    int min_dim = 2;
    bool check_opt_traversal = false;
    for (int i = 0; i < num_tests; ++i)
    {
        int dim1 = rand() % max_dim + min_dim;
        int dim2 = rand() % max_dim + min_dim;
        test_pair_curves(dim1, dim2, check_opt_traversal);
    }
}

void test_dfd_and_traversal(void)
{
    int num_tests = 50;
    int max_dim = 10;
    int min_dim = 2;
    bool check_opt_traversal = true;
    for (int i = 0; i < num_tests; ++i)
    {
        int dim1 = rand() % max_dim + min_dim;
        int dim2 = rand() % max_dim + min_dim;
        test_pair_curves(dim1, dim2, check_opt_traversal);
    }
}

void test_pair_curves(int dim1, int dim2, bool check_opt_traversal)
{
    Curve *c1 = create_curve(dim1);
    Curve *c2 = create_curve(dim2);

    double dist_bf = get_dfd_bruteforce(*c1, 0, *c2, 0);
    double dist_dfd = Metrics::Discrete_Frechet::distance(*c1, *c2);

    if (check_opt_traversal)
    {
        auto lp = std::list<std::tuple<uint32_t, uint32_t>>();
        Metrics::Discrete_Frechet::optimal_traversal(*c1, *c2, lp);
        double dfd_trav = check_traversal(*c1, *c2, lp);
        TEST_ASSERT(dist_bf == dfd_trav);
    }

    TEST_ASSERT(dist_bf == dist_dfd);

    Metrics::Discrete_Frechet::clean();
    delete c1;
    delete c2;
}

Curve *create_curve(int num_points) {
    static int id_counter = 0;
    auto *vec_with_points = new std::vector<Point *>();
    for (int i = 0; i < num_points; i++) {
        double rand_num = rand() % 100 + 1;
        double x_item = rand() % 100 + 1;
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

/* COMPLETE BINARY TREE TESTS */

void test_cbt_post_order_sum(void)
{
    for (int i = 1; i < 1001; i += (rand() % 5 + 1))
    {
        int res_cbt = inner_test_cbt(i);
        int res_act = (i * (i + 1))/2;
        TEST_ASSERT(res_cbt == res_act);
    }
}

int inner_test_cbt(int max_value)
{
    auto *vec = new std::vector<int *>();

    int num_ints = max_value;
    for (int i = 0; i < num_ints; i++) {
        vec->push_back( new int(i+1) );
    }
   
    auto *cbt = new CompleteBinaryTree<int>(vec, vec->size());
    int value = post_order(*cbt);

    for (CBTree_Node node = cbt->get_root(); !cbt->is_leaf(node); ++node) {
        int *item = cbt->get_item(node);
        if (item)
            delete cbt->get_item(node);
    }

    delete cbt;

    for (int i = 0; i < num_ints; i++) {
        TEST_ASSERT(*(vec->at(i)) == (i+1));
    }

    for (auto *v : *vec) { delete v; }
    delete vec;
    return value;
}


int post_order(CompleteBinaryTree<int> &tree)
{
    CBTree_Node root = tree.get_root();
    int *result = internal_post_order(tree, root);
    return *result;
}

int *internal_post_order(CompleteBinaryTree<int> &tree, uint32_t node)
{
    if (tree.is_leaf(node))
        return tree.get_item(node);
    
    int *left_item  = internal_post_order(tree, tree.get_left_child(node));
    int *right_item = nullptr;
    CBTree_Node right_child = tree.get_right_child(node);

    if (!tree.is_empty(right_child))
        right_item = internal_post_order(tree, right_child);

    int *ret_val;
    if (!left_item)
        ret_val = nullptr;
    else if (!right_item)
        ret_val = new int (*left_item);
    else
        ret_val = new int (*left_item + *right_item);

    tree.set_item(node, ret_val == left_item ? nullptr : ret_val );
    return ret_val;
}
