

#ifndef SEARCH_HPP
#define SEARCH_HPP

#include <string>
#include <cstdint>

enum search_type {
    EUCLIDEAN_LSH,
    EUCLIDEAN_HC,
    DISCRETE_FRECHET_LSH,
    CONTINUOUS_FRECHET_LSH
};


class Search {

public:
    Search(enum search_type search_t, std::string &input_path, std::string &query_path, std::string &out_path, double delta=0.5, int k_LSH=5, int k_HC=14, int L=5, int M=10, int probes=2);

    static enum search_type convert_algo_to_type(std::string &algo, std::string &metric) {
        if (algo == "LSH")
        {
            return EUCLIDEAN_LSH;
        }
        else if (algo == "Hypercube")
        {
            return EUCLIDEAN_HC;
        }
        else
        {
            if (metric == "discrete")
                return DISCRETE_FRECHET_LSH;
            else
                return CONTINUOUS_FRECHET_LSH;
        }
    }

};

#endif


