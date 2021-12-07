#include "../src/util/dataset/dataset.hpp"

int main(int argc, char const *argv[]) {

    std::string path = "../data/testdata.txt";
    auto dataset = Dataset(path);

    dataset.print();

    for(auto & curve: *dataset.getData()) {
        curve->filter(2);
        curve->erase_time_axis();
        curve->fit_to_grid(2);

        dataset.print();

        curve->min_max_filter();
        // padding
    }

    dataset.print();
    return 0;
}

