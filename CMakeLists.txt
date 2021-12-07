cmake_minimum_required(VERSION 3.10)
project(project_1)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(lsh tests/main_lsh.cpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/lsh/lsh.cpp src/lsh/lsh.hpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/lsh/hash_functions/amplified_hf.cpp src/lsh/hash_functions/amplified_hf.hpp src/util/utilities.cpp src/util/utilities.hpp src/bruteforce/brute_force_knn.cpp src/bruteforce/brute_force_knn.hpp src/lsh/hashtable/hashtable.hpp src/lsh/hashtable/hashtable.cpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp)
add_executable(cube tests/main_hc.cpp src/hypercube/hypercube.cpp src/hypercube/hypercube.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/util/utilities.cpp src/util/utilities.hpp src/bruteforce/brute_force_knn.cpp src/bruteforce/brute_force_knn.hpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp)
add_executable(cluster tests/main_cluster.cpp src/clustering/cluster.cpp src/clustering/cluster.hpp src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/lsh/hash_functions/amplified_hf.cpp src/lsh/hash_functions/amplified_hf.hpp  src/lsh/lsh.cpp src/lsh/lsh.hpp src/hypercube/hypercube.cpp src/hypercube/hypercube.hpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/utilities.cpp src/util/utilities.hpp src/lsh/hashtable/hashtable.hpp src/lsh/hashtable/hashtable.cpp tests/main_cont_frechet.cpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp)
add_executable(cont_frechet tests/main_cont_frechet.cpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/utilities.cpp src/util/utilities.hpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp src/curve/grid.cpp src/curve/grid.hpp src/curve/curve.cpp src/curve/curve.hpp src/curve/point.hpp)