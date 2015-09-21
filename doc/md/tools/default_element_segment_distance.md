default\_segment\_distance: Euclidean Distance between a Point and a Line Segment
===================

Trajectories are often modelled as sequences of line segments and therefore, a basic
geometric operation is given by calculating the distance between a point and the nearest
point on some segment.

Basically, the operation of finding the distance between a point and a line segment is given by
two steps: 
1.	Finding the nearest point on the line segment and
2.	Calculating the distance

For the first problem, three cases exist: The reference point on the line is either one of the endpoints or a point inside the line segment. In order to calculate this point, an orthogonal
projection from the point to the infinite line is calculated and then, one of the three cases is used to return the actual distance.


Quick Reference
-------------------
	// Class template
	template<class sample_type>
		class default_segment_distance: public element_segment_distance<sample_type,double>	
	// Typical Use
	trajcomp::default_element_segment_distance<std::vector<double>> d;


Sample Code
-------------
The following code sample shows how to instantiate an Euclidean segment distance object and how to calculate the distance between a segment and three points. Note that we are using some modern C++ features, which have to be activated in some compilers.

[default\_element\_segment\_distance.cpp](default_element_segment_distance.cpp)

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

This example first creates two vectors `u` and `v` defining the line segment and then three points:

1.	a = (1,5) will be compared to `u` resulting in a distance of 0.5
2.	b = (3,1.5) will be compared to the segment. The result of the projection is given by (2.25,2.25) for a distance of approximately 1.0607
3.	c = (5,4.5) will be compared to `v`

This sample was contributed by Eike-Jens Hofmann <hoffmannei@cip.ifi.lmu.de>. 


Compilation
------------
In order to compile this sample, libtrajcomp must be installed or the compiler needs to be enabled to find trajcomp.hpp included from the first line. For G++, you can use

	g++ -std=c++11 -o sample default_element_segment_distance.cpp

Note the `-std=c++11` switch, which enables modern C++ support. 

Running the Sample
-----------------
This sample outputs the following text string:

	0.5
	1.06066
	2.5



Concepts and Custom Types
------------------------------

The implementation of the element distance expects a vector-like type. So if you want to use your own types, the following methods must be implemented

* 	`size()` returning the number of elements of your own storage class
*	access operator [] for accessing elements as sample_type
*	double-valued substraction of two objects

You can use the sample at the end of the description of the [Euclidean distance](default_element_distance.html)
for orientation:

