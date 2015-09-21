trajectory: The main class for storing trajectories
===================
When working with spatial data, a very basic operation is to calculate the distance between two points. 
For the trajectory computing library, we decided to represent points as vectors of some user-defined 
type such that you can decide whether these should be double, integer or whatever. For this case, we 
provide the interface `element_distance` and an Euclidean implementation `default_element_distance`.

The Euclidean distance is defined to be the square root of the squares of the distances of the differences of the various dimensions of two vectors.


Quick Reference
-------------------
	// Template
	template <typename ValueType>
	class trajectory : public std::vector<std::vector<ValueType>>

	// Public Members / Types
	unsigned int dimension;
	typedef std::vector<std::vector<ValueType>> Base;
	typedef std::vector<ValueType> ElementType;
	
	// Methods					
	trajectory()
	trajectory(size_t dim)
	void dump()
	void summary()
	load(std::string filename)
	void save(std::string filename)
	};


Details
==============
Basic Idea
--------------
The trajectory class is inherited from a nested vector: The outer dimension contains the
different samples / points of the trajectory while the inner vectors give the individual coordinates
of these points. Subclassing from a nested vector has the consequence, that most algorithms work
on nested vectors even when the class trajectory is not used.

Additionally to the functionality of a nested STL vector, the trajectory class has a public member
dimension, which can be set via a constructor. The method 'push_back' has been redefined to throw
an `std::runtime_error` in case of adding points to a trajectory, which have a different number of elements as compared to the member variable.

Output
-----------
* 	The trajectory class provides a method `dump()` which outputs the contents of a trajectory to
	the standard output and can be handy when debugging programs.
*	The function summary outputs only dimension and size of a trajectory
*	load can be used to load a space-separated value file containing a point on each line
* 	save, on the contrary, can be used to save a space-separated file


Sample Code
-------------
The following code sample shows how to instantiate a Euclidean distance object and how to calculate the distance between two points. Note that we are using some modern C++ features, which have to be activated in some compilers.

[trajectory.cpp](trajectory.cpp)

	#include <trajcomp/trajcomp.hpp>
	#include<iostream>

	using namespace std;

	// specialize our types
	typedef trajcomp::trajectory<double> trajectory;
	
	int main(void)
	{
 		trajectory t(2);
		t.load("./data/prague1.dat");
		t.summary();
		
		trajectory t2;
		t2.load("./data/stair.dat");
		t2.summary();
	
		cout << "d(t2,t2) = "<< trajcomp::discrete_frechet(t2,t2)<< endl;
		cout << "d(t ,t2) = "<< trajcomp::discrete_frechet(t,t2)<< endl;
		cout << "d(t2,t ) = "<< trajcomp::discrete_frechet(t2,t)<< endl;
	
		return 0;	
	}


In this example, first, a trajectory-type based on storing double coordinates is defined to have
the name `trajectory`. Then, a two-dimensional trajectory `t` is created and loaded from the file 
`./data/prague1.dat` contained in the trajcomp release. Then a second two-dimensional trajectory 
`t2` is created and similarly loaded from the file `./data/stair.dat`. Each trajectory ouputs a summary of themselves. Then the discrete Frechet distance of several pairs (t2,t2; t,t2; t2,t) is calculated
and output to the terminal.


Compilation
------------
In order to compile this sample, libtrajcomp must be installed or the compiler needs to be enabled to find trajcomp.hpp included from the first line. For G++, you can use

	g++ -std=c++11 -o sample trajectory.cpp

Note the `-std=c++11` switch, which enables modern C++ support. 

Running the Sample
-----------------
This sample outputs the following text string:

	Dimension = 2
	Elements  = 578
	Dimension = 2
	Elements  = 31
	d(t2,t2) = 0
	d(t ,t2) = 1712.11
	d(t2,t ) = 1712.11

