#include "../src/util/dataset/dataset.hpp"
#include "../src/hypercube/hypercube.hpp"
#include "../src/util/metrics/metrics.hpp"
#include "../src/bruteforce/brute_force_nn.hpp"
#include "../src/util/files/file_reader.hpp"

#include <chrono>

#define GET_DURATION(START, END) (std::chrono::duration_cast<std::chrono::nanoseconds>((END) - (START)).count() * 1e-9)
#define GET_CURR_TIME() (std::chrono::high_resolution_clock::now())


int main(int argc, char const *argv[]) {

    std::string input_path = "../data/nasd_input.csv";
    std::string query_path = "../data/nasd_query_small.csv";
    std::string out_path = "../data/out_euclidean.txt";
    auto dataset = Dataset(input_path);
    auto queries = Dataset(query_path);

    // LSH lsh = LSH(dataset, EUCLIDEAN, 1);
    hypercube hc = hypercube(dataset, Metrics::Euclidean::distance);  // TODO add args
        // std::cout << q->get_size() << std::endl;

    std::cout << "\nRunning LSH with Euclidean distance metric: " << std::endl
         << "- Input path: " << input_path << std::endl
         << "- Query path: " << query_path << std::endl
         << "- Output path: " << out_path << std::endl;

    double avg_lsh_time_taken = 0;
    double abg_brutef_time_taken = 0;
    double maf = 0;

    for(auto query: *queries.getData()) {

        string label = query->get_id();
        // query->erase_time_axis();
        // Benchmark LSH K-NN
        auto start = GET_CURR_TIME();
        auto top_hc = std::multimap<double, std::string>();
        auto q = query->flatten();
        hc.knn(q, top_hc);
        auto end = GET_CURR_TIME();
        double lsh_time_taken = GET_DURATION(start, end);

        // Benchmark Brute-force K-NN
        start = GET_CURR_TIME();
        std::tuple<double, string> top_brutef = {std::numeric_limits<double>::max(), "-"};
        bruteforce_nn(*query, dataset.getData(), Metrics::Discrete_Frechet::distance, &top_brutef);
        end = GET_CURR_TIME();
        auto brutef_time_taken = GET_DURATION(start, end);

        avg_lsh_time_taken += (lsh_time_taken / queries.size());
        abg_brutef_time_taken += (brutef_time_taken / queries.size());

        double candidate_maf = std::get<0>(*top_hc.begin()) / std::get<0>(top_brutef);
        if(std::get<0>(*top_hc.begin()) < std::numeric_limits<double>::max() - 10e5 && maf < candidate_maf)
            maf = candidate_maf;

        std::string name = "Hypercube";
        // std::tuple<double, string> top_brutef = {top_hc.begin()->first, top_hc.begin()->second};
        auto p = std::make_tuple( top_hc.begin()->first, top_hc.begin()->second );
        write_data_to_out_file(label, p, top_brutef, out_path, name);
    }
    write_data_to_out_file(avg_lsh_time_taken, abg_brutef_time_taken, maf, out_path);

    std::cout << "\nGoodbye!" << std::endl;

    return 0;
}

