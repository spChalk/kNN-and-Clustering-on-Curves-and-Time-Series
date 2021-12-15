
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

vector<Point *> *zip_with_time_axis(std::vector<double> &y_values) {

    auto curve = new vector<Point *>();
    double i = 1;
    for(auto & value: y_values) {
        auto new_vector = vector<double>{i++, value};
        curve->push_back(new Point( new_vector ));
    }
    return curve;
}

/*
 *  Reads a line representing a vector in 'd' dimensions and returns
 *  a tuple: (point's ID, vector of doubles).
 */
static Curve *curve_builder(unordered_set<string> &set, const string& input) {

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
    std::vector<double> y_values;
    std::string id = token.c_str();

    // Read the rest coordinates
    while(getline(iss, token, '\t')) {
        if(token == "\r") break;
        y_values.push_back(stod(token, nullptr));
    }

    return new Curve(id, zip_with_time_axis(y_values));
}

static vector<Curve *> *extract_data(ifstream &file) {

    // Initialize data
    auto *curves = new vector<Curve *>();
    auto set = unordered_set<string>();

    // Push every point in the vector
    string line;
    while(getline(file, line)) {
        if(line.empty()) break;
        curves->push_back(curve_builder(set, line));
    }

    return curves;
}

/*
 *  Read the given "filename" and return a vector of points.
 *  Each point consists of a tuple: (point's ID, vector of doubles).
 */
vector<Curve *> *file_reader(const string& filename) {

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
                            std::tuple<double, string> & top_lsh,
                            std::tuple<double, string> & top_brutef,
                            const string& out_path,
                            const string& method) {

    ofstream out;
    out.exceptions(ifstream::badbit);
    try {
        out.open(out_path, std::ios_base::app);

        // Query id
        out << "Query: " << query_id << endl;
        out << "Algorithm: " << method << endl;
        out << "Approximate Nearest neighbor: " << std::get<1>(top_lsh) << endl;
        out << "True Nearest neighbor: " << std::get<1>(top_brutef) << endl;
        out << "distanceApproximate: " << std::get<0>(top_lsh) << endl;
        out << "distanceTrue: " << std::get<0>(top_brutef) << endl << endl;

        out.close();

    } catch (const ofstream::failure &err) {
        ostringstream msg;
        msg << "Exception during the opening of " << out_path << endl;
        throw runtime_error(msg.str());
    }
}

// Writes all info necessary for the output file
void write_data_to_out_file(double lsh_time, double brutef_time, double maf, const string& out_path) {

    ofstream out;
    out.exceptions(ifstream::badbit);
    try {
        out.open(out_path, std::ios_base::app);

        out << "tApproximateAverage: " << lsh_time << endl;
        out << "tTrueAverage: " << brutef_time << endl;
        out << "MAF: " << maf << endl;

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


