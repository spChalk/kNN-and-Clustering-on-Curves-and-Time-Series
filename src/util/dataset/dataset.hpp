
#ifndef PROJECT_1_DATASET_HPP
#define PROJECT_1_DATASET_HPP

#include <vector>
#include <array>
#include <tuple>
#include <cassert>
#include "../../curve/curve.hpp"

// Wrapper for data
class Dataset {

private:
    // Data will consist of a vector consisting of curves
    std::vector<Curve *> *curves;

public:

    Dataset(std::string & data_path);

    ~Dataset();
    std::vector<Curve *> *getData() const;
    uint32_t size();

    // TODO REMOVE IT LATER
    void print();
};

#endif //PROJECT_1_DATASET_HPP
