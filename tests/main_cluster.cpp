#include <cmath>
#include <iostream>
#include "../src/clustering/cluster.hpp"
#include "../src/util/metrics/metrics.hpp"

using std::endl;
using std::get;
using std::ostringstream;
using std::runtime_error;

int main(int argc, char const *argv[]) {

    string input_path, config_path, out_path, method;
    method = "Classic";
    bool verbose = false;

    if(input_handle_cluster(argc, argv, &input_path, &config_path, &out_path, &method, &verbose) == -1)
        exit(EXIT_FAILURE);

    if(input_path.empty())
        input_path = get_path("input file (train dataset)");

    if(config_path.empty())
        config_path = get_path("configuration file");

    // Get output file
    if(out_path.empty())
        out_path = config_output_file();
    else
        remove(out_path.c_str());

    Dataset dataset = Dataset(input_path);

    // method = "Hypercube";
    // method = "LSH";
    auto c = cluster(config_path, dataset, Metrics::euclidean, method);

    c.perform_clustering();
    c.write_results_to_file(out_path, verbose);

    return 0;
}
