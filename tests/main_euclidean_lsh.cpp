#include "../src/util/dataset/dataset.hpp"
#include "../src/lsh/lsh.hpp"
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

    LSH lsh = LSH(dataset, EUCLIDEAN, 1);

    std::cout << "\nRunning LSH with Euclidean distance metric: " << endl
         << "- Input path: " << input_path << endl
         << "- Query path: " << query_path << endl
         << "- Output path: " << out_path << endl;

    double avg_lsh_time_taken = 0;
    double abg_brutef_time_taken = 0;
    double maf = 0;

    for(auto query: *queries.getData()) {

        string label = query->get_id();

        // Benchmark LSH K-NN
        auto start = GET_CURR_TIME();
        std::tuple<double, string> top_lsh = {std::numeric_limits<double>::max(), "-"};
        lsh.nearest_neighbor(query, top_lsh);
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

        double candidate_maf = std::get<0>(top_lsh) / std::get<0>(top_brutef);
        if(std::get<0>(top_lsh) < std::numeric_limits<double>::max() - 10e5 && maf < candidate_maf)
            maf = candidate_maf;

        write_data_to_out_file(label, top_lsh, top_brutef, out_path);
    }
    write_data_to_out_file(avg_lsh_time_taken, abg_brutef_time_taken, maf, out_path);

    std::cout << "\nGoodbye!" << endl;

    return 0;
}