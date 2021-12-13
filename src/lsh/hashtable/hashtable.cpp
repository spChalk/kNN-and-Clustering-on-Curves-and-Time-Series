
#include "hashtable.hpp"
#include <tuple>

using std::get;
using std::tuple;

hashtable::hashtable(uint32_t size, amplified_hf *hf)
        : table(new vector< bucket *>()),
          hashfunction(hf) {
    for(uint32_t i = 0; i < size; i++)
        this->table->push_back(new bucket());
}

hashtable::~hashtable() {

    delete this->hashfunction;
    for(auto b: *table) {
        for(auto pair: *b)
            delete pair;
        b->clear();
        delete b;
    }
    table->clear();
    delete table;
}

void hashtable::insert(FlattenedCurve *new_point) {

    uint32_t id = this->hash(new_point);
    uint32_t index = id % this->table->size();
    (*this->table)[index]->push_back(new tuple<uint32_t, FlattenedCurve *>(id, new_point));
}

uint32_t hashtable::hash(FlattenedCurve *new_point) {
    return this->hashfunction->hash(new_point->get_coordinates());
}

uint32_t hashtable::get_tablesize() {
    return this->table->size();
}

bucket *hashtable::get_bucket(uint32_t index) {
    return (*this->table)[index];
}

