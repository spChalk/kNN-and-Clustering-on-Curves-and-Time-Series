
#include "amplified_hf.hpp"
#include "../../util/utilities.hpp"
#include <vector>

using std::vector;

amplified_hf::amplified_hf(uint32_t hf_num, uint32_t window, uint32_t dim)
: hash_family(new vector<hash_function *>()), random_vars(new vector<int>()) {

    // Initialize random variables and hash tables
    for (uint32_t i = 0; i < hf_num; i++) {
        random_vars->push_back(Distributions::uniform<int>(0, INT32_MAX));
        hash_family->push_back(new hash_function(window, dim));
    }
}

amplified_hf::~amplified_hf() {

    uint32_t size = this->hash_family->size();
    for(uint32_t i = 0; i < size; i++) {
        delete (*this->hash_family)[i];
    }
    delete this->hash_family;
    delete this->random_vars;
}

uint32_t amplified_hf::hash(vector<double> *query){

    // Compute the hash value, by linear combination (R * H)
    uint32_t M = UINT32_MAX - 4;
    uint32_t result = 0;
    for(uint32_t i = 0; i < this->hash_family->size(); i++)
        result += mod((long)(*this->random_vars)[i], M) * mod((long)(*this->hash_family)[i]->hash(query), M);

    return result % M;
}
