
#include <iostream>
#include <string>

#include "../src/clustering/cluster.hpp"
#include "../src/util/utilities.hpp"


int main(int argc, char const *argv[])
{
    std::string input_path, config_path, out_path, update_method, assignment_method;
    bool verbose, silh;

    if(input_handle_cluster(argc, argv, &input_path, &config_path, &out_path, &update_method, &assignment_method, &silh, &verbose) == -1)
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

    std::cout << "Going in " << std::endl;
    cluster(config_path, out_path, dataset, assignment_method, update_method, verbose, silh);
    
    return 0;
}
