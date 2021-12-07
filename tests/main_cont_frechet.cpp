#include "../src/curve/curve.hpp"

int main(int argc, char const *argv[]) {

    std::string path = "../data/testdata.txt";
    auto curve = Curve(path);

    curve.filter(25.0);
    curve.erase_time_axis();
    curve.fit_to_grids(10, 1);

    curve.print();
    return 0;
}

