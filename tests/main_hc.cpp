#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "../src/util/dataset/dataset.hpp"
#include "../src/util/files/file_reader.hpp"
#include "../src/util/utilities.hpp"
#include "../src/bruteforce/brute_force_knn.hpp"
#include "../src/hypercube/hypercube.hpp"
#include "../src/util/metrics/metrics.hpp"

using std::cout;
using std::endl;
using std::ifstream;
using std::ostringstream;
using std::runtime_error;

#define GET_DURATION(START, END) (std::chrono::duration_cast<std::chrono::nanoseconds>((END) - (START)).count() * 1e-9)
#define GET_CURR_TIME() (std::chrono::high_resolution_clock::now())

/* Main function */
int main(int argc, char const *argv[])
{
    // Default parameters
    uint32_t k = 14;  // Projected dimensions
    uint32_t M = 10;  // Max points checked
    uint32_t probes = 2;  // Max vertices checked
    uint32_t N = 1;  // Find N nearest neighbors
    uint32_t R = 10000;  // Find neighbors in radius R

    std::string input_path, query_path, out_path;

    if(input_handle_hypercube(argc, argv, &input_path, &query_path,
                    &out_path, &k, &M, &probes, &N, &R) == -1)
        exit(EXIT_FAILURE);

    if(input_path.empty())
        input_path = get_path("input file (train dataset)");

    // Initialize the dataset by giving the input file's path
    cout << "\nBuilding the Hypercube...\n" << endl;
    Dataset dataset = Dataset(input_path);
    hypercube hc = hypercube(dataset, Metrics::euclidean, k, M, probes, N, R);

    bool repeat = false;
    do
    {
        if(query_path.empty() || repeat)
            query_path = get_path("query file (test dataset)");

        Dataset query = Dataset(query_path);

        if(out_path.empty() || repeat)
            out_path = config_output_file();
        else
            remove(out_path.c_str());

        cout << "\nRunning Hypercube with: " << endl
            << "- Input path: " << input_path << endl
            << "- Query path: " << query_path << endl
            << "- Output path: " << out_path << endl
            << "- Number of Projected Dimensions: " << k << endl
            << "- Max points checked limit: " << M << endl
            << "- Max Hypercube vertices checked limit: " << probes << endl
            << "- Nearest neighbours: " << N << endl
            << "- Range search radius: " << R << endl << endl;

        for(auto pair: *query.getData())
        {
            std::string label = pair->get_id();
            auto test = pair->get_coordinates();

            // Benchmark Hypercube K-NN
            auto start = GET_CURR_TIME();
            auto top_n = std::multimap<double, std::string>();
            hc.knn(pair, top_n);
            auto end = GET_CURR_TIME();
            double time_taken = GET_DURATION(start, end);

            // Benchmark Brute-force K-NN
            start = GET_CURR_TIME();
            auto brutef_top_n = std::multimap<double, std::string>();
            bruteforce_knn(N, *test, dataset.getData(), Metrics::euclidean, &brutef_top_n);
            end = GET_CURR_TIME();
            auto brutef_time_taken = GET_DURATION(start, end);

            // Run Range-query
            auto range_points = std::list<std::tuple<Point *, double>>();
            hc.range_search(pair, range_points);

            write_data_to_out_file(label, top_n, time_taken,
                                  brutef_top_n, brutef_time_taken,
                                  range_points, R, out_path, std::string("CUBE"));
        }
    } while ((repeat = ask_user_to_repeat()));
    
    cout << "\nGoodbye!" << endl;
    return 0;
}
