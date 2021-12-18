#include <iostream>

#include "acutest.h"
#include "../src/clustering/complete_binary_tree.hpp"

int *internal_post_order(CompleteBinaryTree<int> &tree, CBTree_Node node);
int post_order(CompleteBinaryTree<int> &tree);
int inner_test_cbt(int max_value);
void test_sum(void);

TEST_LIST = {
    { "test_sum", test_sum },
    { NULL, NULL }
};

void test_sum(void)
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
    delete cbt;

    for (int i = 0; i < num_ints; i++) {
        TEST_ASSERT(*(vec->at(i)) == (i+1));
    }

    delete vec;
    return value;
}


int post_order(CompleteBinaryTree<int> &tree)
{
    CBTree_Node root = tree.get_root();
    int *result = internal_post_order(tree, root);
    // printf(">> Result is: %d\n\n", result);
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
    return ret_val;
}