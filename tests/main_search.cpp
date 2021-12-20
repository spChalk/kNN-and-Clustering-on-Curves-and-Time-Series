// int input_handle_search(int narg, char const *argvect[], string *inf, string *qf, string *outf, uint32_t *L, uint32_t *delta, uint32_t *k_LSH, uint32_t *k_HC, uint32_t *M, uint32_t *probes, std::string *algo, std::string *metr) {

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
#include "../src/search/search.hpp"

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
int main(int argc, char const *argv[])
{
    // Default parameters
    uint32_t L = 1, k_LSH = 1, k_HC = 1, M = 1, probes = 1;
    double delta = 1.0;

    string input_path, query_path, out_path;
    string algo, metr;

    if(input_handle_search(argc, argv, &input_path, &query_path,
                    &out_path, &L, &delta, &k_LSH, &k_HC, &M, &probes, &algo, &metr) == -1)
        exit(EXIT_FAILURE);

    if(input_path.empty())
        input_path = get_path("input file (train dataset)");

    // Initialize the dataset by giving the input file's path
    enum search_type s_type = Search::convert_algo_to_type(algo, metr);
    auto s = Search(s_type, input_path, query_path, out_path, delta, k_LSH, k_HC, L, M, probes);

    cout << "\nGoodbye!" << endl;
    return 0;
}
