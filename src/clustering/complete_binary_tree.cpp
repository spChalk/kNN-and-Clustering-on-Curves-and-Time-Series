#include <cmath>
#include <cstdint>
#include <assert.h>

#include "complete_binary_tree.hpp"

CompleteBinaryTree::CompleteBinaryTree(void **items, uint32_t total_leaves)
{
    uint32_t height = ceil(log2(total_leaves)) + 1;
    this->max_size = pow(2, height);
    this->array = new void*[this->max_size];
    this->last_level_index = this->max_size / 2 - 1;
    this->actual_size = this->last_level_index + total_leaves;
    for (uint32_t i = 0; i != total_leaves; ++i) {
        this->array[this->last_level_index + i] = items[i];
    }
}

CompleteBinaryTree::~CompleteBinaryTree() {
    delete array;
}

bool CompleteBinaryTree::is_leaf(uint32_t index) {
    return index >= this->last_level_index && index <= actual_size;
}

bool CompleteBinaryTree::is_empty(uint32_t index) {
    return index >= actual_size;
}

void *CompleteBinaryTree::get_item(uint32_t index) {
    assert(index <= actual_size);
    return array[index];
}

void CompleteBinaryTree::set_item(uint32_t index, void *item) {
    assert(index <= actual_size);
    this->array[index] = item;
}

uint32_t CompleteBinaryTree::get_root() {
    return 0;
}

uint32_t CompleteBinaryTree::get_left_child(uint32_t index) {
    return 2*index + 1;
}
uint32_t CompleteBinaryTree::get_right_child(uint32_t index) {
    return 2*index + 2;
}