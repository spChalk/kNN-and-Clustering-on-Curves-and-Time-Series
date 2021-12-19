cmake_minimum_required(VERSION 3.10)
project(project_1)

set(CMAKE_CXX_STANDARD 14)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(lsh tests/main_lsh.cpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/lsh/lsh.cpp src/lsh/lsh.hpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/lsh/hash_functions/amplified_hf.cpp src/lsh/hash_functions/amplified_hf.hpp src/util/utilities.cpp src/util/utilities.hpp src/bruteforce/brute_force_nn.cpp src/bruteforce/brute_force_nn.hpp src/lsh/hashtable/hashtable.hpp src/lsh/hashtable/hashtable.cpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp
        src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/utilities.cpp src/util/utilities.hpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp src/curve/grid.cpp src/curve/grid.hpp src/curve/curve.cpp src/curve/curve.hpp src/curve/point.hpp
        cont_frechet_repo/Fred/include/curve.hpp cont_frechet_repo/Fred/src/curve.cpp
        cont_frechet_repo/Fred/include/point.hpp cont_frechet_repo/Fred/src/point.cpp
        cont_frechet_repo/Fred/include/simplification.hpp cont_frechet_repo/Fred/src/simplification.cpp
        cont_frechet_repo/Fred/include/config.hpp cont_frechet_repo/Fred/src/config.cpp
        cont_frechet_repo/Fred/include/frechet.hpp cont_frechet_repo/Fred/src/frechet.cpp
        src/lsh/lsh.cpp src/lsh/lsh.hpp)

add_executable(cube tests/main_hc.cpp src/hypercube/hypercube.cpp src/hypercube/hypercube.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/util/utilities.cpp src/util/utilities.hpp src/bruteforce/brute_force_nn.cpp src/bruteforce/brute_force_nn.hpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp)

add_executable(test_dfd tests/test_dfd.cpp 
                src/util/metrics/metrics.cpp
                src/util/metrics/metrics.hpp
                src/util/utilities.cpp 
                src/util/utilities.hpp 
                src/curve/grid.cpp 
                src/curve/grid.hpp 
                src/curve/curve.cpp
                src/curve/curve.hpp
                src/curve/point.hpp
                cont_frechet_repo/Fred/include/curve.hpp cont_frechet_repo/Fred/src/curve.cpp
                cont_frechet_repo/Fred/include/point.hpp cont_frechet_repo/Fred/src/point.cpp
                cont_frechet_repo/Fred/include/simplification.hpp cont_frechet_repo/Fred/src/simplification.cpp
                cont_frechet_repo/Fred/include/config.hpp cont_frechet_repo/Fred/src/config.cpp
                cont_frechet_repo/Fred/include/frechet.hpp cont_frechet_repo/Fred/src/frechet.cpp
)

add_executable(cluster tests/main_cluster.cpp 
                src/clustering/cluster.cpp 
                src/clustering/cluster.hpp 
                src/clustering/cluster_template.cpp 
                src/clustering/cluster_template.hpp 
                src/clustering/complete_binary_tree.hpp 
                src/lsh/hashtable/hashtable.hpp 
                src/lsh/hashtable/hashtable.cpp
                src/lsh/hash_functions/hash_function.cpp 
                src/lsh/hash_functions/hash_function.hpp 
                src/lsh/hash_functions/amplified_hf.cpp 
                src/lsh/hash_functions/amplified_hf.hpp  
                src/lsh/lsh.cpp 
                src/lsh/lsh.hpp 
                src/hypercube/hypercube.cpp 
                src/hypercube/hypercube.hpp 
                src/util/dataset/dataset.cpp 
                src/util/dataset/dataset.hpp 
                src/util/files/file_reader.cpp 
                src/util/files/file_reader.hpp 
                src/util/utilities.cpp 
                src/util/utilities.hpp 
                src/util/metrics/metrics.cpp
                src/util/metrics/metrics.hpp
                src/curve/grid.cpp 
                src/curve/grid.hpp 
                src/curve/curve.cpp
                src/curve/curve.hpp
                src/curve/point.hpp
                cont_frechet_repo/Fred/include/curve.hpp cont_frechet_repo/Fred/src/curve.cpp
                cont_frechet_repo/Fred/include/point.hpp cont_frechet_repo/Fred/src/point.cpp
                cont_frechet_repo/Fred/include/simplification.hpp cont_frechet_repo/Fred/src/simplification.cpp
                cont_frechet_repo/Fred/include/config.hpp cont_frechet_repo/Fred/src/config.cpp
                cont_frechet_repo/Fred/include/frechet.hpp cont_frechet_repo/Fred/src/frechet.cpp
        )

# add_custom_target(cluster_test ./cluster -i ../data/nasd_input_small.csv -c ../data/cluster.conf -o ../data/out_cluster_test.txt -update Mean_Vector -assignment Classic)

add_executable(cont_frechet tests/main_cont_frechet.cpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/utilities.cpp src/util/utilities.hpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp src/curve/grid.cpp src/curve/grid.hpp src/curve/curve.cpp src/curve/curve.hpp src/curve/point.hpp
        cont_frechet_repo/Fred/include/curve.hpp cont_frechet_repo/Fred/src/curve.cpp
        cont_frechet_repo/Fred/include/point.hpp cont_frechet_repo/Fred/src/point.cpp
        cont_frechet_repo/Fred/include/simplification.hpp cont_frechet_repo/Fred/src/simplification.cpp
        cont_frechet_repo/Fred/include/config.hpp cont_frechet_repo/Fred/src/config.cpp
        cont_frechet_repo/Fred/include/frechet.hpp cont_frechet_repo/Fred/src/frechet.cpp
        src/lsh/lsh.cpp src/lsh/lsh.hpp
        src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/lsh/hash_functions/amplified_hf.cpp src/lsh/hash_functions/amplified_hf.hpp
        src/bruteforce/brute_force_nn.cpp src/bruteforce/brute_force_nn.hpp
        src/lsh/hashtable/hashtable.hpp src/lsh/hashtable/hashtable.cpp
        )

add_executable(discr_frechet tests/main_discr_frechet.cpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/utilities.cpp src/util/utilities.hpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp src/curve/grid.cpp src/curve/grid.hpp src/curve/curve.cpp src/curve/curve.hpp src/curve/point.hpp
        cont_frechet_repo/Fred/include/curve.hpp cont_frechet_repo/Fred/src/curve.cpp
        cont_frechet_repo/Fred/include/point.hpp cont_frechet_repo/Fred/src/point.cpp
        cont_frechet_repo/Fred/include/simplification.hpp cont_frechet_repo/Fred/src/simplification.cpp
        cont_frechet_repo/Fred/include/config.hpp cont_frechet_repo/Fred/src/config.cpp
        cont_frechet_repo/Fred/include/frechet.hpp cont_frechet_repo/Fred/src/frechet.cpp
        src/lsh/lsh.cpp src/lsh/lsh.hpp
        src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/lsh/hash_functions/amplified_hf.cpp src/lsh/hash_functions/amplified_hf.hpp
        src/bruteforce/brute_force_nn.cpp src/bruteforce/brute_force_nn.hpp
        src/lsh/hashtable/hashtable.hpp src/lsh/hashtable/hashtable.cpp
        )

add_executable(euclidean_lsh tests/main_euclidean.cpp src/util/dataset/dataset.cpp src/util/dataset/dataset.hpp src/util/files/file_reader.cpp src/util/files/file_reader.hpp src/util/utilities.cpp src/util/utilities.hpp src/util/metrics/metrics.cpp src/util/metrics/metrics.hpp src/curve/grid.cpp src/curve/grid.hpp src/curve/curve.cpp src/curve/curve.hpp src/curve/point.hpp
        cont_frechet_repo/Fred/include/curve.hpp cont_frechet_repo/Fred/src/curve.cpp
        cont_frechet_repo/Fred/include/point.hpp cont_frechet_repo/Fred/src/point.cpp
        cont_frechet_repo/Fred/include/simplification.hpp cont_frechet_repo/Fred/src/simplification.cpp
        cont_frechet_repo/Fred/include/config.hpp cont_frechet_repo/Fred/src/config.cpp
        cont_frechet_repo/Fred/include/frechet.hpp cont_frechet_repo/Fred/src/frechet.cpp
        src/lsh/lsh.cpp src/lsh/lsh.hpp
        src/lsh/hash_functions/hash_function.cpp src/lsh/hash_functions/hash_function.hpp src/lsh/hash_functions/amplified_hf.cpp src/lsh/hash_functions/amplified_hf.hpp
        src/bruteforce/brute_force_nn.cpp src/bruteforce/brute_force_nn.hpp
        src/lsh/hashtable/hashtable.hpp src/lsh/hashtable/hashtable.cpp )