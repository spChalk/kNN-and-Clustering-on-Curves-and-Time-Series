#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include "../src/util/dataset/dataset.hpp"
#include "../src/lsh/hash_functions/hash_function.hpp"
#include "../src/lsh/hash_functions/amplified_hf.hpp"
#include "../src/lsh/lsh.hpp"
#include "../src/util/files/file_reader.hpp"
#include "../src/bruteforce/brute_force_nn.hpp"
#include "../src/util/utilities.hpp"
#include "../src/util/metrics/metrics.hpp"

using std::string;
using std::get;
using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ostringstream;
using std::runtime_error;

#define GET_DURATION(START, END) (std::chrono::duration_cast<std::chrono::nanoseconds>((END) - (START)).count() * 1e-9)
#define GET_CURR_TIME() (std::chrono::high_resolution_clock::now())

/* Main function */
int main(int argc, char const *argv[]) {

    // Default parameters
/*    uint32_t num_hf = 4;            // Number of hash functions in each Hash family
    uint32_t num_ht = 5;            // Number of hash tables
    uint32_t num_of_nearest_n = 1;  // Number of nearest neighbours to find
    uint32_t radius = 10000;        // Find neighbours inside a radius

    string input_path, query_path, out_path;

    if(input_handle_lsh(argc, argv, &input_path, &query_path,
                    &out_path, &num_hf, &num_ht, &num_of_nearest_n, &radius) == -1)
        exit(EXIT_FAILURE);

    if(input_path.empty())
        input_path = get_path("input file (train dataset)");

    // Initialize the dataset by giving the input file's path
    cout << "\nBuilding LSH...\n" << endl;
    Dataset dataset = Dataset(input_path);
    // Initialize LSH structure
    LSH lsh = LSH(dataset.getData(), Metrics::euclidean, num_ht, num_hf, num_of_nearest_n, radius);

    bool repeat = false;
    do {
        // Get input file
        if(query_path.empty() || repeat)
            query_path = get_path("query file (test dataset)");

        // Load query data
        Dataset query = Dataset(query_path);

        // Get output file
        if(out_path.empty() || repeat)
            out_path = config_output_file();
        else
            remove(out_path.c_str());

        cout << "\nRunning LSH with: " << endl
             << "- Input path: " << input_path << endl
             << "- Query path: " << query_path << endl
             << "- Output path: " << out_path << endl
             << "- Number of hash functions: " << num_hf << endl
             << "- Number of hash tables: " << num_ht << endl
             << "- Nearest neighbours: " << num_of_nearest_n << endl
             << "- Range search radius: " << radius << endl << endl;

        for(auto pair: *query.getData()) {

            string label = pair->get_id();
            auto test = pair->get_coordinates();

            // Benchmark LSH K-NN
            auto start = GET_CURR_TIME();
            auto top_n = multimap<double, string>();
            lsh.knn(pair, top_n);
            auto end = GET_CURR_TIME();
            double time_taken = GET_DURATION(start, end);

            // Benchmark Brute-force K-NN
            start = GET_CURR_TIME();
            auto brutef_top_n = multimap<double, string>();
            bruteforce_knn(num_of_nearest_n, *test, dataset.getData(), Metrics::euclidean, &brutef_top_n);
            end = GET_CURR_TIME();
            auto brutef_time_taken = GET_DURATION(start, end);

            // Run Range-query
            auto range_points = list<tuple<Point *, double>>();
            lsh.range_search(pair, range_points);

            write_data_to_out_file(label, top_n, time_taken,
                                   brutef_top_n, brutef_time_taken,
                                   range_points, radius, out_path);
        }
    } while ((repeat = ask_user_to_repeat()));

    cout << "\nGoodbye!" << endl;
    return 0;*/
}
