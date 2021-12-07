#include "../src/util/dataset/dataset.hpp"

int main(int argc, char const *argv[]) {

    std::string path = "../data/testdata.txt";
    auto dataset = Dataset(path);

    dataset.print();

    for(auto & curve: *dataset.getData()) {
        curve->filter(2);
        curve->erase_time_axis();
        //TODO: δ must be computed automatically
        curve->fit_to_grid(2);

        dataset.print();

        curve->min_max_filter();
        //TODO: see if padding num is good (put INT16_MAX for testing)
        curve->apply_padding(10);
    }

    dataset.print();
    return 0;
}

