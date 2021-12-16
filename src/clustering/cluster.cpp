


#include "cluster.hpp"
#include "cluster_template.hpp"

#include "../util/dataset/dataset.hpp"
// #include "../lsh/lsh.hpp"
// #include "../hypercube/hypercube.hpp"
#include "../util/utilities.hpp"
#include "../util/files/file_reader.hpp"
#include "../util/metrics/metrics.hpp"

class hypercube;
class lsh;  // TODO : Rm diz -> test


typedef double(*curve_distance_func)(Curve &, Curve &);
typedef double(*flatn_distance_func)(FlattenedCurve &, FlattenedCurve &);

/////////////////////////////////////////////////////////////////////
// Constructor - Destructor(s)

cluster::cluster(std::string &config_path, std::string &out_path, Dataset &dataset, const std::string& assignment, const std::string &update, bool verbose, bool evaluation)
{
    uint32_t _num_of_clusters = 1;
    uint32_t _num_of_ht = 3;
    uint32_t _num_of_hf = 4;
    uint32_t _max_num_of_M_hypercube = 10;
    uint32_t _num_of_hypercube_dims = 3;
    uint32_t _num_of_probes = 2;// Defaults

    // Read the configuration file
    read_cluster_config_file(config_path, &_num_of_clusters, &_num_of_ht, &_num_of_hf,
                             &_max_num_of_M_hypercube, &_num_of_hypercube_dims, &_num_of_probes);

    // Choose the update method
    uint8_t upd_method;
    if (update == "Mean_Vector") {
        upd_method = 0;
    }
    else if (update == "Mean_Frechet") {
        upd_method = 1;
        // this->lsh_ds = new LSH(this->dataset.getData(), this->dist_func, _num_of_ht, _num_of_hf);  // TODO : hmm
    }
    else {
        ostringstream msg;
        msg << "Error, Update Method should be one of the following:\n> Mean_Vector\n> Mean_Frechet\n." << endl;
        throw runtime_error(msg.str());
    }

    uint8_t assignment_method;
    // Choose the execution method
    if (assignment == "Classic") {
        assignment_method = 0;
    }
    else if (assignment == "LSH") {
        assignment_method = 1;
        // this->lsh_ds = new LSH(this->dataset.getData(), this->dist_func, _num_of_ht, _num_of_hf);
    }
    else if (assignment == "LSH_Frechet") {
        assignment_method = 1;
        // this->lsh_ds = new LSH(this->dataset.getData(), this->dist_func, _num_of_ht, _num_of_hf);  // TODO: hmm
    }
    else if (assignment == "Hypercube") {
        assignment_method = 2;
        // this->hc_ds = new hypercube(this->dataset, this->dist_func, _num_of_hypercube_dims, _max_num_of_M_hypercube, _num_of_probes);
    }
    else {
        ostringstream msg;
        msg << "Error, Assignment Method should be one of the following:\n> Classic\n> LSH\n> Hypercube\n> LSH_Frechet\n." << endl;
        throw runtime_error(msg.str());
    }

    if (upd_method == 0)  // Vectors <-> FlattenedCurves
    {
        std::vector<FlattenedCurve *> *data_to_cluster = dataset.flatten_data();
        if (assignment_method == 0)  // Classic Lloyd's
        {
            auto cluster_obj = internal_cluster<FlattenedCurve, flatn_distance_func>(_num_of_clusters, data_to_cluster, &Metrics::Euclidean::distance, assignment_method, upd_method, nullptr, nullptr);
            cluster_obj.perform_clustering();
            cluster_obj.write_results_to_file(out_path, verbose, evaluation);
        }
        else if (assignment_method == 1)  // LSH - Initialize for "FlattenedCurve" dataset
        {

        }
        else  // Hypercube for "FlattenedCurve" dataset
        {

        }
        // delete data_to_cluster ?
    }
    else  // Curves
    {
        std::vector<Curve *> *data_to_cluster = dataset.getData();
        if (assignment_method == 0)  // Classic Lloyd's
        {
            auto cluster_obj = internal_cluster<Curve, curve_distance_func>(_num_of_clusters, data_to_cluster, &Metrics::Discrete_Frechet::distance, assignment_method, upd_method, nullptr, nullptr);
            cluster_obj.perform_clustering();
            cluster_obj.write_results_to_file(out_path, verbose, evaluation);
        }
        else  // LSH - Initialize for "Curve" dataset
        {

        }
    }
}
