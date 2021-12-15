
#ifndef PROJECT_1_FILE_READER_HPP
#define PROJECT_1_FILE_READER_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include "../../curve/curve.hpp"

/*
 *  Read the given "filename" and return a vector of points.
 *  Each point consists of a tuple: (point's ID, vector of doubles).
 */
std::vector<Curve *> *file_reader(const std::string& filename);

// Writes all info necessary for the output file
void write_data_to_out_file(const std::string& query_id,
                            std::tuple<double, std::string>&,
                            std::tuple<double, std::string>&,
                            const std::string&, const std::string& method="LSH_Frechet_Continuous");

void write_data_to_out_file(double lsh_time, double brutef_time, double maf, const std::string& out_path);

void read_cluster_config_file(const std::string& filename, uint32_t *noc, uint32_t *noht,
                              uint32_t *nohf, uint32_t *mMh, uint32_t *nhd, uint32_t *nop);

#endif //PROJECT_1_FILE_READER_HPP
