#include "trajcomp/trajcomp.hpp"

int main() {
    std::vector<double> u {1.0d, 1.0d};
	std::vector<double> v {3.0d, 3.0d};
    std::vector<double> a {1.0d, 0.5d};
    std::vector<double> b {3.0d, 1.5d};
    std::vector<double> c {5.0d, 4.5d};

    trajcomp::default_segment_distance<std::vector<double>> d;
	    
    // Should be 0.5
    std::cout << d(u, v, a) << std::endl;
    // Should be sqrt(1.125) = 1.0607
    std::cout << d(u, v, b) << std::endl;
    // Should be 2.5
    std::cout << d(u, v, c) << std::endl;
    return 0;
}
