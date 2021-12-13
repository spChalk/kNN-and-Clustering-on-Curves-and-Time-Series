#include "../src/util/dataset/dataset.hpp"
#include "../src/lsh/lsh.hpp"
#include "../src/util/metrics/metrics.hpp"
#include "../src/util/utilities.hpp"

int main(int argc, char const *argv[]) {

    std::string input_path = "../data/nasd_input_small.csv";
    std::string query_path = "../data/nasd_query.csv";
    auto dataset = Dataset(input_path);
    auto queries = Dataset(query_path);

    uint32_t grid_interval = estimate_grid_interval(dataset, queries);
    uint32_t max_curve_length = compute_max_curve_length(dataset, queries);

    for(auto & curve: *dataset.getData()) {
        curve->filter(2);
        curve->erase_time_axis();
        curve->fit_to_grid(grid_interval);
        curve->min_max_filter();
        //TODO: see if padding num is good (put INT16_MAX for testing)
        curve->apply_padding(max_curve_length);
    }

    auto lsh_data = dataset.flatten_data();
    LSH lsh = LSH(lsh_data, Metrics::continuous_frechet_distance, 1);
    for (auto &d: *lsh_data)
        delete d;
    delete lsh_data;

    return 0;
}

