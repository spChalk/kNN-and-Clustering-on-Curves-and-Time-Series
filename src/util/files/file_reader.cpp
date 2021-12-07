
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <iomanip>
#include "file_reader.hpp"

using std::endl;
using std::tuple;
using std::string;
using std::vector;
using std::multimap;
using std::list;
using std::unordered_set;
using std::unordered_multimap;
using std::make_tuple;
using std::get;
using std::stod;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::ostringstream;
using std::runtime_error;

/*
 *  Reads a line representing a vector in 'd' dimensions and returns
 *  a tuple: (point's ID, vector of doubles).
 */
static Point *vector_builder(unordered_set<string> &set, const string& input) {

    istringstream iss(input);
    string token;

    // Read the first token (ID)
    getline(iss, token, '\t');

    /* If the ID is already present in the duplicates set,
     * throw an error and exit. */
    if(set.find(token) != set.end()) {
        ostringstream msg;
        msg << "Error, input has duplicate IDs: " << token.c_str() << endl;
        throw runtime_error(msg.str());
    } else {
        // Else insert the ID in the set
        set.insert(token);
    }

    // Build the tuple and insert its unique ID
    std::vector<double> coord;
    std::string id = token.c_str();

    // Read the rest coordinates
    while(getline(iss, token, '\t')) {
        if(token == "\r") break;
        coord.push_back(stod(token, nullptr));
    }

    return new Point(id, coord);
}

static vector<Point *> *extract_data(ifstream &file) {

    // Initialize data
    auto *data = new vector<Point *>();
    auto set = unordered_set<string>();

    // Push every point in the vector
    string line;
    while(getline(file, line)) {
        if(line.empty()) break;
        data->push_back(vector_builder(set, line));
    }

    return data;
}

/*
 *  Read the given "filename" and return a vector of points.
 *  Each point consists of a tuple: (point's ID, vector of doubles).
 */
vector<Point *> *file_reader(const string& filename) {

    ifstream f(filename.c_str());
    if (!f.good())
    {
        ostringstream msg;
        msg << "File in path: \"" << filename << "\" doesn't exist." << endl;
        throw runtime_error(msg.str());
        exit(1);
    }

    ifstream file;
    file.exceptions(ifstream::badbit);
    try {

        // Open the file
        file.open(filename);

        // Initialize data
        auto *data = extract_data(file);

        // Close the file and return the data
        file.close();
        return data;

    } catch (const ifstream::failure &err) {
        ostringstream msg;
        msg << "Exception during the opening of " << filename << endl;
        throw runtime_error(msg.str());
    }
}

// Writes all info necessary for the output file
void write_data_to_out_file(const string& query_id,
                            multimap<double, string> & lsh_data,
                            double method_time_taken,
                            multimap<double, string> & brutef_data,
                            double brutef_time_taken,
                            list<tuple<Point *, double>> & range_points,
                            int radius,
                            const string& out_path,
                            const string& method) {

    ofstream out;
    out.exceptions(ifstream::badbit);
    try {
        out.open(out_path, std::ios_base::app);

        // Query id
        out << "Query: " << query_id << endl << endl;

        // Compare LSH to Hypercube
        uint32_t i = 1;
        for(auto it_1 = lsh_data.cbegin(), end_1 = lsh_data.cend(),
                    it_2 = brutef_data.cbegin(), end_2 = brutef_data.cend();
            it_1 != end_1 || it_2 != end_2;) {

            if(it_1 != end_1) {
                out << "Nearest neighbor-" << i++ << ": " << it_1->second << endl;
                out << "distance" << method << ": " << it_1->first << endl;
            }
            if(it_2 != end_2) {
                if(it_1 == end_1) {
                    out << "Nearest neighbor-" << i++ << ": " << it_2->second << endl;
                }
                out << "distanceTrue: " << it_2->first << endl << endl;
                ++it_2;
            }
            if(it_1 != end_1) ++it_1;
        }

        // Print time intervals
        out << "t" << method << ": " << method_time_taken << std::setprecision(5) << endl;
        out << "tTrue: " << brutef_time_taken << std::setprecision(5) << endl;

        // Print all points inside the provided radius
        out << radius << "-near neighbors:" << endl;
        for(const auto& point: range_points)
            out << (*(get<0>(point))).get_id() << endl;
        out << endl;

        out.close();

    } catch (const ofstream::failure &err) {
        ostringstream msg;
        msg << "Exception during the opening of " << out_path << endl;
        throw runtime_error(msg.str());
    }
}

void read_cluster_config_file(const string& filename, uint32_t *noc, uint32_t *noht,
                              uint32_t *nohf, uint32_t *mMh, uint32_t *nhd, uint32_t *nop) {
    ifstream f(filename.c_str());
    if (!f.good())
    {
        ostringstream msg;
        msg << "File in path: \"" << filename << "\" doesn't exist." << endl;
        throw runtime_error(msg.str());
        exit(1);
    }

    ifstream file;
    file.exceptions(ifstream::badbit);
    try {

        // Open the file
        file.open(filename);

        string line;
        while(getline(file, line)) {

            if(line.empty()) break;
            istringstream iss(line);
            string token;

            // Read the first token (ID)
            getline(iss, token, ' ');
            auto label = string(token);

            // Read the rest coordinates
            while(getline(iss, token, ' ')) {
                if(token == "\r") break;
                if(label == "number_of_clusters:")
                    *noc = strtol(token.c_str(), nullptr, 10);
                else if(label == "number_of_vector_hash_tables:")
                    *noht = strtol(token.c_str(), nullptr, 10);
                else if(label == "number_of_vector_hash_functions:")
                    *nohf = strtol(token.c_str(), nullptr, 10);
                else if(label == "max_number_M_hypercube:")
                    *mMh = strtol(token.c_str(), nullptr, 10);
                else if(label == "number_of_hypercube_dimensions:")
                    *nhd = strtol(token.c_str(), nullptr, 10);
                else if(label == "number_of_probes:")
                    *nop = strtol(token.c_str(), nullptr, 10);
            }
        }

        // Close the file and return the data
        file.close();

    } catch (const ifstream::failure &err) {
        ostringstream msg;
        msg << "Exception during the opening of " << filename << endl;
        throw runtime_error(msg.str());
    }
}


