
#ifndef PROJECT_1_HASHTABLE_HPP
#define PROJECT_1_HASHTABLE_HPP

#include <vector>
#include <tuple>
#include "../../util/dataset/dataset.hpp"
#include "../hash_functions/amplified_hf.hpp"

using std::vector;
using std::tuple;
using std::string;

typedef vector< tuple<uint32_t , Point *> *> bucket;

/*
 * Custom constant size hash table implemented with separate chaining
 */
class hashtable {

private:
    vector< bucket *> *table;
    amplified_hf *hashfunction;

public:
    hashtable(uint32_t size, amplified_hf *hf);
    virtual ~hashtable();
    void insert(Point *);
    uint32_t hash(Point *);
    uint32_t get_tablesize();
    bucket *get_bucket(uint32_t);
};

#endif //PROJECT_1_HASHTABLE_HPP
