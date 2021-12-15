
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>
#include <functional>
#include "hypercube.hpp"

#include "../util/utilities.hpp"  // window calculation

static std::mt19937 rng;  // Mersenne twister random number generator
static auto rnd = std::bind(std::uniform_int_distribution<uint32_t>(0,1),std::ref(rng));  // Random number drawn uniformly (p(0) = p(1) = 0.5) from {0, 1}

// Set limits related to query operations.
void hypercube::set_limits(uint32_t newN, uint32_t newR, uint32_t newM, uint32_t new_probes) {
    this->M = newM; this->probes = new_probes; this->N = newN; this->R = newR;
}
// TODO : Test init dim
// Constructor - creates a Hypercube data structure and inserts every record from given dataset.
hypercube::hypercube(Dataset &dataset, distance_f dist_func, uint32_t k, uint32_t M, uint32_t probes, uint32_t N, uint32_t R)
: init_dim(dataset.getData()->at(0)->get_data_dimensions()), k(k), dist_func(dist_func)
{
    this->set_limits(N, R, M, probes);

    uint32_t num_buckets = pow(2, k);  // 2^(d') buckets
    this->hypercube_ht = new vector<vector<FlattenedCurve *>>(num_buckets);
    this->bucket_to_vertex_index = new std::unordered_map<uint64_t , int>();
    this->hash_family = new vector<hash_function *>(k);

    auto data = dataset.flatten_data();
    uint32_t window = estimate_window_size(data, this->dist_func);

    for (uint32_t i = 0; i < k; i++) {  // Create d' LSH-admitting hash functions
        (*(this->hash_family))[i] = new hash_function(window, init_dim);
    }
    for (auto point_ptr: *data)
        this->insert_point(point_ptr);  // Insert dataset records
    
    // delete data;  // TODO: Hopefully this will not invoke item destructors
}

hypercube::~hypercube() {
    delete this->hypercube_ht;
    delete this->bucket_to_vertex_index;
    for (uint32_t i = 0; i < this->k; i++) {
        delete (*(this->hash_family))[i];
    }
    delete this->hash_family;
}

void hypercube::insert_point(FlattenedCurve *point_ptr) {
    uint32_t bucket = convert_bucket_to_vertex(point_ptr);
    (*(this->hypercube_ht))[bucket].push_back(point_ptr);
}

// Project a bucket given by h-functions to a Hypercube vertex.
uint32_t hypercube::convert_bucket_to_vertex(FlattenedCurve *datapoint)
{
    vector<double> *coordinates = datapoint->get_coordinates();
    uint32_t vertex = 0x0000;

    // For every h function
    for (uint32_t i = 0; i < this->hash_family->size(); ++i)
    {
        vertex <<= 1;  // Shift left since each h_i is stored in 1 bit of the vertex
        uint32_t bucket = (*(this->hash_family))[i]->hash(coordinates);  // Get h_i(p)

        // key=(i, bucket), unique key so we can remember which pair of (h_i, bucket) we have already projected
        uint64_t key = (((uint64_t) i) << 32) + (uint64_t) bucket;

        // Check if we already mapped (h_i, h(p)) to {0, 1} in the past
        uint32_t mapping;
        auto found = bucket_to_vertex_index->find(key);
        if (found == bucket_to_vertex_index->end())  // Not found
        {
            mapping = rnd(rng);  // Generate uniformly 0 or 1
            bucket_to_vertex_index->insert( {key, mapping} );  // Remember this choice
        } else {
            mapping = found->second;  // Already calculated in the past - get its mapping
        }
        // Construct the vertex by adding the mapping for h_i
        vertex += mapping;  // mapping can only change the last bit of the vertex
    }
    return vertex;
}

// Find the vertex of a point in the Hypercube
uint32_t hypercube::search_for_vertex(FlattenedCurve *p) {
    return convert_bucket_to_vertex(p);
}

// Pick the next vertex to explore
uint32_t hypercube::select_next_vertex(bool & visited_all, uint32_t current_vertex, uint32_t max_vertices, std::queue<uint32_t> &next_vertices, std::unordered_set<uint32_t> & vertex_cache)
{
    vertex_cache.insert(current_vertex);  // Remember that we have explored current vertex
    if (vertex_cache.size() == max_vertices) {
        visited_all = true;
        return -1;  // We explored every vertex of the Hypercube, thus stop.
    }

    // Find every vertex with HD=1 wrt the current vertex, and insert them in a FIFO queue.
    // We do that by flipping exactly 1 bit of the current vertex binary representation in each loop.
    for (uint32_t i = 0; i < max_vertices; ++i)
    {
      uint32_t neighb = current_vertex ^ ((uint32_t) 0x0001 << i);  // XOR just 1 bit <-> HD=1

      if (vertex_cache.find(neighb) == vertex_cache.end())
        next_vertices.push(neighb);  // Check if we already explored "neighb" vertex
    }

    uint32_t next_vertex = next_vertices.front();  // Pick the next vertex
    next_vertices.pop();
    return next_vertex;
}

// Update the list storing the results of range queries.
static void maintain_result_container(std::list<std::tuple<FlattenedCurve *, double>> & top_n, uint32_t radius, double dist, std::string & label, FlattenedCurve *p)
{
    if (dist < radius)
        top_n.emplace_back(p, dist);
}

// Update the multimap storing the K results of a K-NN query.
static void maintain_result_container(std::multimap<double, std::string> & top_n, uint32_t k, double dist, std::string & label, FlattenedCurve *p)
{
    if (top_n.size() == k) {
        if (dist < top_n.rbegin()->first)
            top_n.erase(--top_n.end());
        else
            return;
    }
    top_n.insert({dist, label});
}

// Perform a K-NN query.
void hypercube::knn(FlattenedCurve *query, std::multimap<double, std::string> & top_n) {
  if(this->N == 0) return;
  this->perform_query(query, this->N, this->dist_func, top_n);
}

// Perform a range query.
void hypercube::range_search(FlattenedCurve *query, std::list<std::tuple<FlattenedCurve *, double>>& top_n) {
  this->perform_query(query, this->R, this->dist_func, top_n);
}

// Perform a query (either K-NN or Range)
// C denotes the type of the Container storing the Results of the Query
// T denotes the metric used for bounds check (type of k and R parameters)
template <typename C, typename T>
void hypercube::perform_query(FlattenedCurve *query, T metric, distance_f dist_func, C & top_n)
{
    // Keep explored vertices
    auto vertex_cache = std::unordered_set<uint32_t>();
    // Keep vertices with priority (FIFO) to be explored next
    auto next_vertices = std::queue<uint32_t>();
    
    uint32_t points_checked = 0;
    uint32_t vertices_checked = 0;

    // Find the vertex where query is mapped
    uint32_t vertex = this->search_for_vertex(query);
    bool visited_all_v = false;

    while (!visited_all_v)
    {
        auto bucket_v = &(*(this->hypercube_ht))[vertex];

        // For every (train) point in the current bucket
        for(auto it : *bucket_v)
        {
            std::string label = it->get_id();
            // Calculate distance and check if it belongs in the result
            double dist = dist_func(*query, *it);
            maintain_result_container(top_n, metric, dist, label, it);
            if (++points_checked >= this->M)
                break;
        }

        if (++vertices_checked == this->probes)
            break;  // Check if the stop conditions are met

        vertex = select_next_vertex(visited_all_v, vertex, this->k, next_vertices, vertex_cache);
    }
}
