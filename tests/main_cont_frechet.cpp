#include "../src/util/dataset/dataset.hpp"

int main(int argc, char const *argv[]) {

    std::string path = "../data/testdata.txt";
    auto dataset = Dataset(path);

    /*dataset.filter(25.0);
    dataset.erase_time_axis();
    dataset.fit_to_grids(10, 1);
*/
    dataset.print();
    return 0;
}

