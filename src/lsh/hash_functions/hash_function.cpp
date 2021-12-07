
#include "hash_function.hpp"
#include "../../util/utilities.hpp"
#include <cmath>

hash_function::hash_function(uint32_t window, uint32_t dim):
noise(Distributions::uniform<double>(0, window)),
normal_vector(new vector<double>()) {

    for(uint32_t i = 0; i < dim; i++)
        normal_vector->push_back(Distributions::normal<double>());

    divide_vector_by_scalar<double>(this->normal_vector, window);
    this->noise /= (double)window;
}

int hash_function::hash(vector<double> *query) {

    return floor(dot_product<double>(query, this->normal_vector) + this->noise);
}

hash_function::~hash_function() {
    delete normal_vector;
}
