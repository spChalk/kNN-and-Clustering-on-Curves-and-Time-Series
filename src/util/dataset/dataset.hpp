
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
    std::vector<FlattenedCurve *> *flatten_data();
    std::vector<FlattenedCurve *> *erase_time_and_flatten_data();
    uint32_t size();

    void print();
};

#endif //PROJECT_1_DATASET_HPP
