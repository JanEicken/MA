#include<trajcomp/trajcomp.hpp>
#include<vector>

using namespace std;
using namespace trajcomp;

int main() {
    std::vector<double> p {1.0d, 1.0d};
    std::vector<double> q {3.0d, 3.0d};

    trajcomp::default_element_distance<std::vector<double>> d;

cout << "The distance between "
     << tools::make_string(p) 
     << " and "
     << tools::make_string(q) 
     << " is "
     << d(p,q) << endl; 

    return 0;
}
