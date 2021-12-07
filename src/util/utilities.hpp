
#ifndef PROJECT_1_UTILITIES_HPP
#define PROJECT_1_UTILITIES_HPP

#include <iostream>
#include <random>
#include <vector>
#include <sstream>
#include "dataset/dataset.hpp"

using std::endl;
using std::string;
using std::vector;
using std::ostringstream;
using std::runtime_error;

class Distributions {

public:

    // Uniform distribution for integers and real numbers
    template<typename T>
    static T uniform(long from, long to) {
        static std::random_device                  rand_dev;
        static std::mt19937                        mt_generator(rand_dev());
        return std::is_integral<T>::value ?
                std::uniform_int_distribution<int>(from, to)(mt_generator)
                : std::uniform_real_distribution<double>(from, to)(mt_generator);
    }

    // Normal distribution
    template<typename T>
    static T normal() {
        static std::random_device         rand_dev;
        static std::default_random_engine r_generator(rand_dev());
        std::normal_distribution<T>     distribution(0, 1);
        return distribution(r_generator);
    }
};

template<typename T>
void divide_vector_by_scalar(vector<T> *v, uint32_t scalar) {
    for(uint32_t i = 0; i < v->size(); i++)
        v->data()[i] = v->data()[i] / scalar;
}

template<typename T>
void multiply_vector_by_scalar(vector<T> *v, uint32_t scalar) {
    for(uint32_t i = 0; i < v->size(); i++)
        v->data()[i] = v->data()[i] * scalar;
}

template<typename T>
void add_scalar_to_vector(vector<T> *v, double scalar) {
    for(uint32_t i = 0; i < v->size(); i++)
        v->data()[i] = v->data()[i] + scalar;
}

template<typename T>
void add_vectors(vector<T> *a, vector<T> *b, vector<T> *result) {

    if(a->size() != b->size() || result->size() != a->size()) {
        ostringstream msg;
        msg << "Exception in 'add_vectors'. Vectors do not have the same size." << endl;
        throw runtime_error(msg.str());
    }

    for(uint32_t i = 0; i < a->size(); i++)
        result->data()[i] = a->data()[i] + b->data()[i];
}

template<typename T>
int vectors_are_equal(vector<T> *a, vector<T> *b) {
    if(a->size() != b->size()) {
        ostringstream msg;
        msg << "Exception in 'add_vectors'. Vectors do not have the same size." << endl;
        throw runtime_error(msg.str());
    }
    for(uint32_t i = 0; i < a->size(); i++) {
        double diff = a->data()[i] - b->data()[i];
        if (diff > 1e-5 || diff < -1e-5)
            return 0;
    }
    return 1;
}

template<typename T>
double dot_product(vector<T> *a, vector<T> *b) {

    if(a->size() != b->size()) {
        ostringstream msg;
        msg << "Exception in 'dot_product'. Vectors do not have the same size. Size A = "
            << a->size() << ", Size B = " << b->size() << endl;
        throw runtime_error(msg.str());
    }

    double prod = 0;
    for(uint32_t i = 0; i < a->size(); i++)
        prod = prod + a->data()[i] * b->data()[i];
    return prod;
}

/*
 * Estimates the value of the window size in LSH and Hypercube implementations.
 * 'w' should be approx. 1 - 10 times the average distance between all points.
 * Because this app handles big data, the average distance is being computed between
 * a subset (1%) of data.
 */
uint32_t estimate_window_size(vector<Point *> *data, double(*distance)(Point&, Point&));

// Modulo operation between two numbers
uint32_t mod(long a, uint32_t b);

// Gets a file specified by the user and returns it
string get_path(const string& mode);

// Gets a file specified by the user, and returns it (removes the existing one if exists)
string config_output_file();

// Receives a string and decides whether the string is an integer or not.
bool is_uint (char *str);

bool ask_user_to_repeat();

// Fully handles the LSH input provided to the app
int input_handle_lsh(int narg, char const *argvect[], string *inf, string *qf, string *outf, uint32_t *k, uint32_t *L, uint32_t *NN, uint32_t *radius);

// Fully handles the Hypercube input provided to the app
int input_handle_hypercube(int narg, char const *argvect[], string *inf, string *qf, string *outf, uint32_t *k, uint32_t *M, uint32_t *probes, uint32_t *NN, uint32_t *radius);

// Fully handles the Cluster input provided to the app
int input_handle_cluster(int narg, char const *argvect[], string *inf, string *cf, string *outf, string *method, bool *verbose);

// Fully copies a point and returns the newly allocated one
Point *copy_point(Point *p);

void apply_floor_to_vector(vector<double> *v);

#endif //PROJECT_1_UTILITIES_HPP
