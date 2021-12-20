
#include "search.hpp"
#include "../util/dataset/dataset.hpp"
#include "../lsh/lsh.hpp"
#include "../hypercube/hypercube.hpp"
#include "../util/utilities.hpp"
#include "../util/files/file_reader.hpp"
#include "../util/metrics/metrics.hpp"
#include "../bruteforce/brute_force_nn.hpp"

#include <chrono>

#define GET_DURATION(START, END) (std::chrono::duration_cast<std::chrono::nanoseconds>((END) - (START)).count() * 1e-9)
#define GET_CURR_TIME() (std::chrono::high_resolution_clock::now())


Search::Search(enum search_type search_t, std::string &input_path, std::string &query_path, std::string &out_path, double delta, int k_LSH, int k_HC, int L, int M, int probes)
{
    if (search_t == EUCLIDEAN_LSH)
    {   
        auto dataset = Dataset(input_path);
        auto queries = Dataset(query_path);
        LSH lsh = LSH(dataset, CONTINUOUS_FRECHET, L, k_LSH);  // Must be L=1 !

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
    }
    else if (search_t == EUCLIDEAN_HC)
    {
        auto dataset = Dataset(input_path);
        auto queries = Dataset(query_path);

        hypercube hc = hypercube(dataset, Metrics::Euclidean::distance, k_HC, M, probes);  // TODO add args

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
    }
    else if (search_t == DISCRETE_FRECHET_LSH)
    {
        auto dataset = Dataset(input_path);
        auto queries = Dataset(query_path);

        LSH lsh = LSH(dataset, DISCRETE_FRECHET, L, k_LSH);

        std::cout << "\nRunning LSH with Discrete distance metric: " << endl
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

            write_data_to_out_file(label, top_lsh, top_brutef, out_path, "LSH_Discrete_Frechet");
        }
        write_data_to_out_file(avg_lsh_time_taken, abg_brutef_time_taken, maf, out_path);

        std::cout << "\nGoodbye!" << endl;
    }
    else
    {
        auto dataset = Dataset(input_path);
        auto queries = Dataset(query_path);

        LSH lsh = LSH(dataset, CONTINUOUS_FRECHET, 1, k_LSH);  // Must be L=1 !

        std::cout << "\nRunning LSH with Discrete distance metric: " << endl
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

            write_data_to_out_file(label, top_lsh, top_brutef, out_path, "LSH_Continuous_Frechet");
        }
        write_data_to_out_file(avg_lsh_time_taken, abg_brutef_time_taken, maf, out_path);

        std::cout << "\nGoodbye!" << endl;
    }
}