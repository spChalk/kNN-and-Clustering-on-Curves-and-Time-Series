#ifndef COMPLETE_BINARY_TREE_TMPL_H
#define COMPLETE_BINARY_TREE_TMPL_H

#include <cmath>
#include <cstdint>
#include <vector>
#include <assert.h>

typedef uint32_t CBTree_Node;

template <typename ItemType>
class CompleteBinaryTree
{
private:
    std::vector<ItemType *> *array;
    uint32_t actual_size;       // Max total items in the CBT
    uint32_t last_level_index;  // Index of the first leaf

public:
    CompleteBinaryTree(std::vector<ItemType *> *leaves, uint32_t total_leaves)
    {
        uint32_t height = ceil(log2(total_leaves)) + 1;
        uint32_t max_size = pow(2, height);
        this->last_level_index = max_size / 2 - 1;
        this->actual_size = this->last_level_index + total_leaves;
        this->array = new std::vector<ItemType *>(this->actual_size, nullptr);
        for (uint32_t i = 0; i != total_leaves; ++i) {
            this->array->at(this->last_level_index + i) = leaves->at(i);
        }
    }

    ~CompleteBinaryTree() {
        delete this->array;
    }
    
    bool is_leaf(CBTree_Node index) {
        return index >= this->last_level_index;
    }

    bool is_internal(CBTree_Node index) {
        return index < this->last_level_index;
    }

    bool is_empty(CBTree_Node index) {
        return index >= actual_size;
    }

    ItemType *get_item(CBTree_Node index) {
        if (index < actual_size)
            return this->array->at(index);
        return nullptr;
    }

    void set_item(CBTree_Node index, ItemType *item) {
        assert(index <= actual_size);
        this->array->at(index) = item;
    }

    CBTree_Node get_root() {
        return 0;
    }

    CBTree_Node get_left_child(CBTree_Node index) {
        return 2*index + 1;
    }

    CBTree_Node get_right_child(CBTree_Node index) {
        return 2*index + 2;
    }
};

#endif