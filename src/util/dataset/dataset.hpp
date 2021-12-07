
#ifndef PROJECT_1_DATASET_HPP
#define PROJECT_1_DATASET_HPP

#include <vector>
#include <array>
#include <tuple>
#include <cassert>

/*
 * A point on space is a tuple consisting of:
 *  1. Point's label
 *  2. A vector of doubles (actual coordinates)
 */
class Point {

    private:
        std::string id;
        std::vector<double> *coordinates;

    public:

        Point(std::string &_id, std::vector<double> &_coord)
        : id(_id), coordinates(new std::vector<double>(_coord)) {}

        static Point *copy_point(Point &p) {
            return new Point(p.id, (*p.coordinates));
        }

        ~Point() {
            delete this->coordinates;
        }

        std::string &get_id() {
            return this->id;
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

};

// Wrapper for data
class Dataset {

protected:
    // Data will consist of a vector consisting of points
    std::vector<Point *> *data;

public:

    Dataset(std::string & data_path);

    virtual ~Dataset();
    std::vector<Point *> *getData() const;
    uint32_t get_data_dimensions();
    uint32_t size();
};

#endif //PROJECT_1_DATASET_HPP
