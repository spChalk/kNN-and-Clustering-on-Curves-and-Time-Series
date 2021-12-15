
#include <cstdint>

class CompleteBinaryTree
{
private:
    void **array;
    uint32_t max_size;
    uint32_t actual_size;
    uint32_t last_level_index;

public:
    bool is_leaf(uint32_t index);
    bool is_empty(uint32_t index);
    void *get_item(uint32_t index);
    void set_item(uint32_t, void *item);
    uint32_t get_root();
    uint32_t get_left_child(uint32_t index);
    uint32_t get_right_child(uint32_t index);

    CompleteBinaryTree(void **leaves, uint32_t total_leaves);
    ~CompleteBinaryTree();
};

