#include<trajcomp/trajcomp.hpp>
#include<vector>
#include<iostream>
using namespace std;
int main(void)
{
  std::vector<double> A={1,2,3,4,5};
  cout << "make_string: " << trajcomp::tools::make_string(A) << endl;
  cout << "make_string (csv): " << trajcomp::tools::make_string<std::vector<double>,','>(A) << endl;
  return 0;
}
