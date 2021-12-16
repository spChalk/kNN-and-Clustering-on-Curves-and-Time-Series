#include "../src/util/dataset/dataset.hpp"
#include "../src/lsh/lsh.hpp"
#include "../src/util/metrics/metrics.hpp"
#include "../src/bruteforce/brute_force_nn.hpp"

int main(int argc, char const *argv[]) {

    std::string input_path = "../data/nasd_input.csv";
    std::string query_path = "../data/nasd_query.csv";
    std::string out_path = "../data/out_discr_frechet.txt";
    auto dataset = Dataset(input_path);
    auto queries = Dataset(query_path);

    LSH lsh = LSH(dataset, DISCRETE_FRECHET, 2);

    std::cout << "\nRunning LSH with Discrete Frechet distance metric: " << endl
         << "- Input path: " << input_path << endl
         << "- Query path: " << query_path << endl
         << "- Output path: " << out_path << endl;

    std::tuple<double, string> result = {std::numeric_limits<double>::max(), "-"};
    auto query = (*queries.getData())[0];
    lsh.nearest_neighbor(query, result);

    std::cout << "\nGoodbye!" << endl;

    return 0;
}

