
#ifndef PROJECT_2_CLUSTER_HPP
#define PROJECT_2_CLUSTER_HPP

#include <string>
#include "../util/dataset/dataset.hpp"

class cluster {

public:
    cluster(std::string &config_path, 
            std::string &out_path,
            Dataset &dataset, 
            const std::string &assignment, 
            const std::string &update,
            bool verbose,
            bool evaluation);
};

#endif //PROJECT_2_CLUSTER_HPP
