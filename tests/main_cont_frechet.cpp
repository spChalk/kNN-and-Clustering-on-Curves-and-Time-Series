#include "../src/util/dataset/dataset.hpp"
#include "../src/lsh/lsh.hpp"
#include "../src/util/metrics/metrics.hpp"

int main(int argc, char const *argv[]) {

    std::string path = "../data/nasd_input_small.csv";
    auto dataset = Dataset(path);

    for(auto & curve: *dataset.getData()) {
        curve->filter(2);
        curve->erase_time_axis();
        //TODO: δ must be computed automatically
        curve->fit_to_grid(2);
        curve->min_max_filter();
        //TODO: see if padding num is good (put INT16_MAX for testing)
        curve->apply_padding(500);
    }

    auto lsh_data = dataset.flatten_data();
    LSH lsh = LSH(lsh_data, Metrics::continuous_frechet_distance, 1);
    for (auto &d: *lsh_data)
        delete d;
    delete lsh_data;

    return 0;
}

