#include<trajcomp/trajcomp.hpp>
#include<vector>
using namespace std;
using namespace trajcomp::tools;

typedef vector< vector <double>> matrix2d;
typedef vector< vector <vector<double>>> matrix3d;


int main(void)
{
	matrix2d A;
	// Create a 3x3 matrix
	matrix_resize(A,3,3);
	matrix3d B;
	// Create a 3x3x2 matrix
	matrix_resize(B,3,3,2);
	return 0;
}
