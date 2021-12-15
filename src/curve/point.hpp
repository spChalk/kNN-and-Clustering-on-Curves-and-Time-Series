
#ifndef PROJECT_1_POINT_HPP
#define PROJECT_1_POINT_HPP

#include <vector>
#include <cstdint>
#include <cassert>
#include <iostream>

/*
 * A point on space is a tuple consisting of:
 *  1. Point's label
 *  2. A vector of doubles (actual coordinates)
 */
class Point {

private:
    std::vector<double> *coordinates;

public:

    Point(std::vector<double> &_coord)
            : coordinates(new std::vector<double>(_coord)) {}

    ~Point() {
        delete this->coordinates;
    }

    double get_coordinate(uint32_t index) {
        assert (index < get_dimensions());
        return this->coordinates->at(index);
    }

    std::vector<double> *get_coordinates() {
        return this->coordinates;
    }

    uint32_t get_dimensions() {
        return coordinates->size();
    }

    // TODO REMOVE IT WHEN READY
    void print() {
        std::cout << "( ";
        for(auto i: *coordinates) {
            std::cout << i << " ";
        }
        std::cout << ")\n";
    }

};

#endif //PROJECT_1_POINT_HPP
