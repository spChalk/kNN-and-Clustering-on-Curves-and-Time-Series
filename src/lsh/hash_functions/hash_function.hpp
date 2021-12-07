
#ifndef PROJECT_1_HASH_FUNCTION_HPP
#define PROJECT_1_HASH_FUNCTION_HPP

#include <vector>
#include <cstdint>

using std::vector;

// LSH's hash function
class hash_function {

private:
    double noise;
    vector<double> *normal_vector;

public:
    hash_function(uint32_t window, uint32_t dim);
    virtual ~hash_function();

    /*
     * h(p) = floor( (p * v + t) / w ) =
     *        floor( (p*v) / w + t/w ) =
     *        floor( p * (v/w) + (t/w) ) =
     *        floor( p * a + b )
     */
    int hash(vector<double> *query);
};

#endif //PROJECT_1_HASH_FUNCTION_HPP
