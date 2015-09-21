matrix_resize: Using nested vectors to store matrices
===================

Sometimes, the trajectory computing library needs 2D and 3D matrices. For compatibility with
iterator-based algorithms, we decided to use nested STL vectors as matrices. However,
creating regular matrices (of a fixed dimension) is tedious. To simplify this, we provide
two function templates, helping allocating matrices.


Quick Reference
-------------------
	template<class datatype> 
	void matrix_resize(std::vector<std::vector< datatype> > &m, int r, int c)

	template<class datatype> 
	void matrix_resize(std::vector<std::vector<std::vector< datatype > > > &m, int r, int c, int d)

Sample Code
-------------
The following code sample shows how to use the matrix class:

[matrix_resize.cpp](matrix_resize.cpp)

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

This example creates two matrices of given dimensions by looping over the different container.

Compilation
------------
In order to compile this sample, libtrajcomp must be installed or the compiler needs to be enabled to find trajcomp.hpp included from the first line. For G++, you can use

	g++ -std=c++11 -o sample matrix_resize.cpp

Note the `-std=c++11` switch, which enables modern C++ support. 

Running the Sample
-----------------
The sample does not create useful output.

Concepts and Custom Types
------------------------------

The implementations expect vector-like types which must provide:

* 	`resize()` method for allocating space
*	access operator [] for accessing elements as data_type



