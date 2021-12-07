
#ifndef PROJECT_1_AMPLIFIED_HF_HPP
#define PROJECT_1_AMPLIFIED_HF_HPP

#include <vector>
#include <cstdint>
#include "hash_function.hpp"

using std::vector;

// LSH's amplified hash function
class amplified_hf {

private:
    // Family of hash functions
    vector<hash_function *> *hash_family;
    // Natural random variables
    vector<int> *random_vars;

public:
    amplified_hf(uint32_t hf_num, uint32_t window, uint32_t dim);
    virtual ~amplified_hf();
    uint32_t hash(vector<double> *query);
};

#endif //PROJECT_1_AMPLIFIED_HF_HPP
