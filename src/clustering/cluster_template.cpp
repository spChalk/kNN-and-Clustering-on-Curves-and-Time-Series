
#include "cluster_template.hpp"

using std::tuple;
using std::stringstream;
using std::vector;
using std::map;
using std::unordered_map;
using std::unordered_multimap;
using std::set;
using std::get;

/////////////////////////////////////////////////////////////////////
// Constructor - Destructor(s)

template <typename ItemType, typename DistFunc>
internal_cluster<ItemType, DistFunc>::
internal_cluster(uint32_t _num_clusters, std::vector<ItemType *> *_dataset, DistFunc _dist_func, uint8_t asgn_method, uint8_t updt_method, LSH *_lsh_ds, hypercube *_hc_ds)
:
dataset(_dataset),
assign_method(asgn_method),
update_method(updt_method),
num_of_clusters(_num_clusters),
dist_func(_dist_func),
centroids(new set<ItemType *>()),
lsh_ds(_lsh_ds),
hc_ds(_hc_ds),
final_assign(new unordered_map<ItemType *, vector<ItemType *> *>())
{ }

// Fully delete the app's centroids
template <typename ItemType, typename DistFunc>
void 
internal_cluster<ItemType, DistFunc>::
delete_centroids(set< ItemType *> *centroid_set) {
    if (centroid_set == nullptr)
        return;
    for(auto centroid: *centroid_set) {
        delete centroid;
    }
    delete centroid_set;
}

// Clear the app's assignment (does not delete assignment's pointer)
template <typename ItemType, typename DistFunc>
void 
internal_cluster<ItemType, DistFunc>::
delete_assignment() {
    for(auto & item: *this->final_assign) {
        item.second->clear();
        delete item.second;
    }
    this->final_assign->clear();
}

template <typename ItemType, typename DistFunc>
internal_cluster<ItemType, DistFunc>::
~internal_cluster() {
    this->delete_centroids(this->prev_centroids);
    this->delete_centroids(this->centroids);
    this->delete_assignment();
    delete this->final_assign;
    delete this->lsh_ds;
    delete this->hc_ds;
}
/////////////////////////////////////////////////////////////////////
// ASSIGNMENTS

// Exact Lloyd's assign
// Assigns every datapoint to its exact nearest centroid
template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
assign_exact_lloyds()
{
    for(auto datapoint: *(this->dataset)) {

        // Initialize minimum distance to a large value
        double min_dist = DBL_MAX;
        ItemType *nearest_centroid = nullptr;
        // For every centroid
        for(auto centroid: *(this->centroids)) {
            // Compute the respective distances and choose the current datapoint's closest centroid
            double dist = this->dist_func(*datapoint, *centroid);
            if(dist < min_dist) {
                min_dist = dist;
                nearest_centroid = centroid;
            }
        }
        // Update assignments structure
        assign_point_to_centroid(datapoint, nearest_centroid, this->final_assign);
    }
}

template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
assign_point_to_centroid(ItemType *pt, ItemType *ctr, std::unordered_map<ItemType *, std::vector<ItemType *> *> *final_assign)
{
    auto centroid_entry = final_assign->find(ctr);

    // if the centroid has 0 assigned points, initialize its vector
    if (centroid_entry == final_assign->end()) {
        auto assigned_points = new std::vector<ItemType *>();
        assigned_points->push_back(pt);
        final_assign->insert( {ctr, assigned_points} );
    }
    else
        centroid_entry->second->push_back(pt);
}

// Get the minimum distance between 2 centroids
template <typename ItemType, typename DistFunc>
double
internal_cluster<ItemType, DistFunc>::
get_min_dist_between_2_centroids()
{
    // Initialize minimum distance to a large value
    double min_dist = DBL_MAX, dist;
    // Loop efficiently through all centroids and fetch the minimum distance
    for (auto it1 = this->centroids->begin(); it1 != this->centroids->end(); ++it1) {
        for (auto it2 = it1; it2 != this->centroids->end(); ++it2)
        {
            if (*it1 != *it2)
            {
                dist = this->dist_func(*(*it1), *(*it2));
                if (dist < min_dist)
                    min_dist = dist;
            }
        }
    }
    return min_dist;
}


// C denotes the underlying `container` to perform range search queries
template <typename ItemType, typename DistFunc>
template <typename C>
void
internal_cluster<ItemType, DistFunc>::
reverse_assignment(C cont)
{
    // Insert every point from the dataset vector in the `unassigned_points` set
    auto unassigned_points = std::unordered_set<ItemType *>(dataset->begin(), dataset->end());

    // Get initial radius
    double radius = this->get_min_dist_between_2_centroids();
    bool little_ball_updates = false;
    uint32_t MIN_UNASS_PTS = 87;

    // Perform reverse assignment, doubling the radius in each iteration
    // until the vast majority of points are assigned or most balls get no updates
    while (unassigned_points.size() >= MIN_UNASS_PTS && !little_ball_updates)
    {
        // Keep track of all the points we assigned for a *fixed* radius iteration
        auto curr_rad_assign =  std::unordered_map<ItemType *, std::tuple<ItemType *, double>>(); 
        // For every centroid
        for (auto centroid : *(this->centroids))
        {
            // Perform a range search with `centroid` as its base
            auto result_list = std::list<std::tuple<ItemType *, double>>();
            cont->range_search(centroid, result_list);

            for (auto &result_tuple : result_list)  // For every point in the range
            {
                ItemType *datapoint = std::get<0>(result_tuple);

                // Ff datapoint is unassigned in the *final* assignment (i.e. wasn't assigned earlier)
                if (unassigned_points.find(datapoint) != unassigned_points.end())
                {
                    double dist_from_curr_centroid = std::get<1>(result_tuple);
                    auto rad_centroid_asgn = curr_rad_assign.find(datapoint);
                    // Check if it was already assigned to another centroid in the *same* fixed radius
                    if (rad_centroid_asgn == curr_rad_assign.end())
                    {
                        curr_rad_assign.insert( {datapoint, std::make_tuple(centroid, dist_from_curr_centroid)} );
                    }
                    else  // ItemType already assigned for current fixed radius
                    {
                        // Compare distance to previous assigned centroid vs distance of current centroid found
                        if (std::get<1>(rad_centroid_asgn->second) > dist_from_curr_centroid)
                        {   // Update and map "datapoint" to the new centroid, if distance is smaller
                            rad_centroid_asgn->second = std::make_tuple(centroid, dist_from_curr_centroid);
                        }
                    }
                }
            }
        }

        // Keep track of the balls (defined by `centroids`) updated for this radius
        std::list<ItemType *> balls_updated;

        // Convert *current radius assignments* made, to *final assignments*
        for (auto &entry : curr_rad_assign)
        {
            auto point_to_assign = entry.first;
            auto centroid_assigned = std::get<0>(entry.second);
            assign_point_to_centroid(point_to_assign, centroid_assigned, this->final_assign);
            unassigned_points.erase(point_to_assign);  // Remove ItemType from `unassigned ItemTypes` set
            balls_updated.push_front(centroid_assigned);
        }

        // If most balls (more than half) don't change, end reverse assignment
        if (balls_updated.size() < num_of_clusters / 2)
            little_ball_updates = true;

        // Double the radius for the next iteration
        radius *= 2;
    }

    // Assign the rest of the points via brute force
    for (auto unass_point : unassigned_points)
    {
        double min_dist = DBL_MAX;
        ItemType *nearest_centroid;

        for(auto centroid: *centroids)  // Compare distances to every centroid
        {
            double dist = this->dist_func(*unass_point, *centroid);
            if(dist < min_dist) {
                min_dist = dist;
                nearest_centroid = centroid;
            }
        }
        assign_point_to_centroid(unass_point, nearest_centroid, this->final_assign);
    }
}

template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
assign_LSH() {
    this->reverse_assignment(this->lsh_ds);
}

template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
assign_Hypercube() {
    this->reverse_assignment(this->hc_ds);
}

/////////////////////////////////////////////////////////////////////
// Initialization

// Get the minimum distance between a given point and a centroid
template <typename ItemType, typename DistFunc>
double
internal_cluster<ItemType, DistFunc>::
min_distance_from_centroids(ItemType *curr_point,
                            std::unordered_multimap<string, std::tuple<string, double>>& history)
{
    // Initialize minimum distance to a large value
    double min_dist = DBL_MAX;
    // For every centroid
    for(auto centroid: *centroids) {
        double dist = -1;
        string curr_point_label = curr_point->get_id();
        string centroid_label = centroid->get_id();

        // Find the computed distances in history, for current point
        auto range = history.equal_range(curr_point_label);
        for(auto it = range.first; it != range.second; it++) {
            // If the distance between the current point and the current
            //      centroid has been computed, fetch it
            if(get<0>(it->second) == centroid_label) {
                dist = get<1>(it->second);
                break;
            }
        }
        // Else, compute the distance and save it to history
        if(dist < 0) {
            dist = this->dist_func(*curr_point, *centroid);
            history.insert({curr_point_label, {centroid_label, dist}});
        }

        if(dist < min_dist)
            min_dist = dist;
    }
    return min_dist;
}

// Centroids initialization
// Implementation of Initialization++ algorithm
template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
initialize_centroids()
{
    auto data = this->dataset;
    // Pick a random point for initial centroid
    int index = Distributions::uniform<int>(0, (long)data->size());
    this->centroids->insert( ItemType(*((*data)[index])) );  // TODO : check this
                                                             // copy-constructor
     // Write down all the distances computed
     // Map point -> (centroid, distance)
    auto history = std::unordered_multimap<string, tuple<string, double>>();
    // Repeat until there are as many centroids as requested
    while(this->centroids->size() < this->num_of_clusters) {

        unordered_map< ItemType *, double > distances;
        // For all points in the dataset
        for(auto _point: *data) {
            // If current point is not a centroid
            if(this->centroids->find(_point) == this->centroids->end()) {
                // Compute its distance from the nearest centroid
                distances.insert( {_point, min_distance_from_centroids(_point, history)} );
            }
        }

        normalize_distances(distances);
        double denominator = compute_sum_of_squared_distances(distances);

        // Compute the probability of an item to be chosen as the next centroid
        unordered_map< ItemType *, double > probabilities;
        for(auto & dist: distances)
            probabilities.insert( {dist.first, ((dist.second * dist.second) / denominator) * 100} );

        // Pick a number uniformly between 0 and 100
        auto uniform_index = Distributions::uniform<double>(0, 100);

        // Iterate through the map until the cumulative probability of the
        // visited elements is greater than the selected number
        double cum_prob = 0.0;
        for(auto & prob: probabilities) {
            cum_prob += prob.second;
            if(uniform_index < cum_prob) {
                centroids->insert(ItemType(*prob.first));  // TODO : Check this 
                break;              // previously was: copy_point(prob.first)
            }                       // so we probably need a copy-constructor? 
        }
    }
}


/////////////////////////////////////////////////////////////////////
// UPDATE

// Find maximum distance in a map consisting of distances
template <typename ItemType, typename DistFunc>
double
internal_cluster<ItemType, DistFunc>::
max_distance(unordered_map<ItemType *, double> &distances) {
    double max_dist = 0.0;
    for(auto & dist: distances) {
        if(dist.second > max_dist)
            max_dist = dist.second;
    }
    return max_dist;
}

// Normalize a map of distances
template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
normalize_distances(unordered_map<ItemType *, double> &distances) {
    // Fetch and divide every element by the max distance
    double max_dist = max_distance(distances);
    for(auto & dist: distances)
        dist.second /= max_dist;
}

template <typename ItemType, typename DistFunc>
double
internal_cluster<ItemType, DistFunc>::
compute_sum_of_squared_distances(unordered_map<ItemType *, double> &distances) {
    double sumsq = 0.0;
    for(auto & dist: distances)
        sumsq += (dist.second * dist.second);
    return sumsq;
}

// Update function
// Calculate the mean of all vector points per cluster and update the centroids
template <typename ItemType, typename DistFunc>
double
internal_cluster<ItemType, DistFunc>::
update_vector()
{
    auto new_centroids = new set<ItemType *>();
    double max_distance_between_centroids = -1.0;
    uint32_t counter = 0;
    // For each centroid
    for(auto & centroid_family: *this->final_assign) {

        auto centroid = centroid_family.first;
        auto points = centroid_family.second;

        // Initialize an array of zeros
        auto result = new vector<double>(centroid->dimensions, 0.0);
        // Compute the mean of all points in the current centroid
        for(auto & _point: *points) {
            vector<double> temp = *(_point->get_coordinates());
            divide_vector_by_scalar<double>(&temp, points->size());
            add_vectors<double>(result, &temp, result);
        }
        // Create the new, updated centroid  // TODO : Test this
        auto new_centroid = new ItemType(string("centroid_").append(std::to_string(++counter)), *result);
        new_centroids->insert(new_centroid);

        // Get the maximum value from the amount of changes of each centroid
        double temp_max_dist = this->dist_func(*centroid, *new_centroid);
        if(temp_max_dist > max_distance_between_centroids)
            max_distance_between_centroids = temp_max_dist;
    }
    // Delete old centroids
    this->delete_centroids(prev_centroids);
    this->prev_centroids = this->centroids;
    this->centroids = new_centroids;

    return max_distance_between_centroids;
}

// Update function
// Calculate the mean of all vector points per cluster and update the centroids
template <typename ItemType, typename DistFunc>
double
internal_cluster<ItemType, DistFunc>::
update_curve()
{
    return 42.0;
}

/////////////////////////////////////////////////////////////////////
// Runs the clustering algorithm and computes its execution time
template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
perform_clustering()
{
    #define GET_DURATION(START, END) (std::chrono::duration_cast<std::chrono::nanoseconds>((END) - (START)).count() * 1e-9)
    #define GET_CURR_TIME() (std::chrono::high_resolution_clock::now())

    auto assign_func = this->assign_method == 0 ? &internal_cluster<ItemType, DistFunc>::assign_exact_lloyds:
                        this->assign_method == 1 ? &internal_cluster<ItemType, DistFunc>::assign_LSH:
                        &internal_cluster<ItemType, DistFunc>::assign_Hypercube;

    auto update_func = this->update_method == 0 ? &internal_cluster<ItemType, DistFunc>::update_vector
                                                : &internal_cluster<ItemType, DistFunc>::update_curve;
    auto start = GET_CURR_TIME();

    this->internal_run(assign_func, update_func);
    
    auto end = GET_CURR_TIME();
    this->time_taken = GET_DURATION(start, end);
}

// Main running function for clustering
template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
internal_run(__CLUSTER_TMPL_MODULE_assign_func assign_f, __CLUSTER_TMPL_MODULE_update_func update_f)
{
    // Initialize the centroids
    this->initialize_centroids();
    // While the maximum value from the amount of changes of each centroid
    // is more than a specified threshold, repeat the
    // update and assignment processes
    do {
        this->delete_assignment();
        ((*this).*assign_f)();
    }
    while(((*this).*update_f)() > 10);
}

/////////////////////////////////////////////////////////////////////
// EVALUATE - SILHOUETTE

template <typename ItemType, typename DistFunc>
ItemType *
internal_cluster<ItemType, DistFunc>::
get_2nd_closest_centroid(ItemType *current_centroid, ItemType *query)
{
    double max_dist = DBL_MAX;
    ItemType *second_closest = nullptr;
    
    for (auto &cluster : *(this->final_assign))
    {
        auto centroid = std::get<0>(cluster);
        if (centroid != current_centroid)
        {
            double dist = this->dist_func(*query, *centroid);
            if (dist < max_dist)
            {
                max_dist = dist;
                second_closest = centroid;
            }
        }
    }
    return second_closest;
}

// Evaluate Clustering based on the Silhouette metric
template <typename ItemType, typename DistFunc>
double
internal_cluster<ItemType, DistFunc>::
evaluate(std::list<double> &result_per_cluster)
{
    double overall = 0.0;  // Overall Silhouette

    for (auto &cluster1 : *(this->final_assign))  // For every cluster
    {
        double s = 0.0;  // Cluster's Silh. value
        auto cluster_centroid = cluster1.first;
        auto assigned_points = cluster1.second;
        for (auto &point1: *assigned_points)  // For every point in the cluster
        {
            double a = 0;  // Calculate avg distance between points in the same cluster
            for (auto &point2: *assigned_points)
                a += this->dist_func(*point1, *point2);
            a /= assigned_points->size();

            // Get second closest centroid/cluster
            ItemType *second_closest_ctr = this->get_2nd_closest_centroid(cluster_centroid, point1);
            auto found = this->final_assign->find(second_closest_ctr);
            if (found == this->final_assign->end()) {
                continue;
            }

            auto &cluster2 = *found;
            auto assigned_points2 = std::get<1>(cluster2);
            double b = 0;  // Calculate avg distance with points in the 2nd closest cluster
            for (auto point3 : *assigned_points2)
                b += this->dist_func(*point1, *point3);
            b /= assigned_points2->size();

            double si = (b - a) / ( a > b ? a : b );  // point's Silh. value
            s += si;
            overall += si;
        }
        s /= assigned_points->size();
        result_per_cluster.push_back(s);
    }

    overall /= this->dataset->size();
    return overall;
}

static void write_vec_to_file(std::ofstream &out, std::vector<double> &vec)
{
    uint32_t dims = vec.size();
    for (uint32_t i=0; i < dims; ++i)
    {
        if (i == 0) out << "[ ";
        out << vec[i];
        if (i != dims-1)
            out << ", ";
        else
            out << " ]";
    }
}


template <typename ItemType, typename DistFunc>
void
internal_cluster<ItemType, DistFunc>::
write_results_to_file(const std::string & out_path, bool verbose, bool evalution_on)
{
    using std::ofstream;
    using std::ostringstream;
    using std::runtime_error;

    ofstream out;
    out.exceptions(std::ifstream::badbit);
    try {
        out.open(out_path, std::ios_base::app);

        // TODO : fix diz
        // // Algorithm
        // if (method == 0)
        //     out << "Algorithm: Lloyds" << endl;
        // else if (method == 1)
        //     out << "Algorithm: Range Search LSH" << endl;
        // else if (method == 2)
        //     out << "Algorithm: Range Search Hypercube" << endl;

        uint32_t index = 0;
        for (auto &cluster : *(this->final_assign))
        {
            auto cluster_centroid = std::get<0>(cluster);
            auto assigned_points = std::get<1>(cluster);
            out << "CLUSTER-" << ++index << " {size: " << assigned_points->size() << ", centroid: ";

            auto centroid_coords = cluster_centroid->get_coordinates();
            write_vec_to_file(out, *centroid_coords);
            out << "}" << std::endl;
        }

        out << "clustering_time: " << time_taken << std::setprecision(5) << " sec" << std::endl;

        if (evalution_on)
        {
            auto result_per_cluster = std::list<double>();
            double overall_sil = this->evaluate(result_per_cluster);

            result_per_cluster.push_back(overall_sil);
            out << "Silhouette: ";

            std::vector<double> res_per_cluster_vec { std::make_move_iterator(std::begin(result_per_cluster)), 
              std::make_move_iterator(std::end(result_per_cluster)) };

            write_vec_to_file(out, res_per_cluster_vec);
            out << std::endl;
        }

        if (verbose)
        {
            out << std::endl;
            uint32_t index = 0;
            for (auto &cluster : *(this->final_assign))
            {
                auto cluster_centroid = std::get<0>(cluster);
                auto assigned_points = std::get<1>(cluster);
                out << "CLUSTER-" << ++index << " {centroid: ";
                auto centroid_coords = cluster_centroid->get_coordinates();
                write_vec_to_file(out, *centroid_coords);

                for (auto p : *assigned_points)
                {
                    auto p_label = p->get_id();
                    auto p_coords = p->get_coordinates();

                    out << ", " << p_label << ": ";
                    write_vec_to_file(out, *p_coords);
                }
                out << "}" << std::endl;
            }
        }

        out.close();

    } catch (const ofstream::failure &err) {
        ostringstream msg;
        msg << "Exception during the opening of " << out_path << endl;
        throw runtime_error(msg.str());
    }
}