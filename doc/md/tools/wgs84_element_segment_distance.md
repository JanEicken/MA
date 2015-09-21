wgs84\_segment\_distance: WGS84 Distance Between a Point and a Line Segment
===================

Trajectories are often modelled as sequences of line segments and therefore, a basic
geometric operation is given by calculating the distance between a point and the nearest
point on some segment. In many applications, however, segments are defined in terms of GPS coordinates for which the shortest connection between two points is not a straight line.

Though there might be more direct ways to deal with this problem, we treat this as a minimization
problem on the set of distances: Find the minimal distance between a point on the geodesic line
forming the `line segment` on the surface of the earth and a fixed point p.

Therefore, we use an iterative optimization method, which avoids calculating derivatives (except one very coarse approximation to this, implicitly). We observe that the minimum distance lies somewhere inbewteen the first and the last point of the geodesic and that only one minimum can exist. Therefore, it is sufficient to do the following, where L denotes the geodesic line between the points and
L(t) denotes the point, when this line is parametrized with a time t from [0,1], such that L(0) is the first point and L(1) is the last point.

Given a Geodesic Line between U and V,

1.	Set t1=0, t4=1 and calculate distances d1=d(Q,L(t1)) and d4=d(Q,L(t4))
2.	Unless the size of the interval (d4-d1) is small enough
	1.	Put two inner points t2,t3 uniformly inbetween t1 and t4
	2. 	Calculate the point on the geodesic for t2 and t3	
	3.	Calculate the distances to Q for the inner points, e.g., d2 = d(Q,L(t2)); d3 = d(Q,L(t3))
	4. 	If d1 is minimal among all distances or d2<d3, refine to the left setting t3 = t4 (and d3 = d4). Otherwise refine to the right setting t1 = t2 (and d1 = d2).
3.	Return the smaller value of d2, d3

This algorithm has proven effective and efficient (though there is room for improvement). However, we decided to take this way in order to 

*	have access to the nearest point on the segment afterwards, which keeps stored in the distance functional object
* 	avoid local numerical differentiation 

Note, howver, calculating the difference of d2 and d3 can be seen as a very rough approximation to the derivative of the distance function in the middle of the interval [t1,t4].


Quick Reference
-------------------
	// Class template
	template<class sample_type>
		class wgs84_segment_distance: public element_segment_distance<sample_type,double>	
	// Typical Use
	trajcomp::wgs84_element_segment_distance<std::vector<double>> dseg;
	dseg.getLastNearestPoint()


Sample Code
-------------
The following code sample shows how to instantiate a WGS84 distance object as well
as a WGS84 point line distance object and how to use these to calculate relations between
three major german cities. 
Note that we are using some modern C++ features, which have to be activated in some compilers.

[wgs84\_distances.cpp](wgs84_distances.cpp)

	#include <trajcomp/trajcomp_geo.hpp>
	using namespace trajcomp;

	typedef std::vector<double> coord;

	int main(void)
	{
	
	    coord munich	=	{48.135125,11.581981};
	    coord nuremberg = 	{49.455555,11.078611};
	    coord cologne   =   {50.937531, 6.960279};
	    wgs84_element_distance<coord> d;
	    wgs84_segment_distance<coord> dseg;
	   
	   cout << "Munich => Nuremberg: "
	        << d(munich, nuremberg) << "m" <<endl;
	
	   cout << "Nuremberg => Cologne: "
	        << d(nuremberg, cologne) << "m" <<endl;

	   cout << "Munich => Cologne: "
	        << d(munich,cologne) << "m" <<endl;

	   cout << 	"Distance of Nuremberg from the "<<endl << 
				"shortest path on earth" << endl <<
			    "(geodesic) between Munich and Cologne: " <<
			    dseg(munich,cologne,nuremberg) << endl;
	
	   cout << "Point on Munich=>Cologne: "
	        <<  tools::make_string(dseg.getLastNearestPoint())
	        <<endl;
		    
	
    return 0;
}


Compilation
------------
In order to compile this sample, libtrajcomp must be installed or the compiler needs to be enabled to find trajcomp.hpp included from the first line. For G++, you can use

	g++ -std=c++11 -lGeographic -o wgs84 wgs84_distances.cpp 


The `-std=c++11` switch enables modern C++ support, the `-lGeographic` switch links with the geodesic libary in use.

Running the Sample
-----------------
This sample outputs the following text string:

	Munich => Nuremberg: 151424m
	Nuremberg => Cologne: 337043m
	Munich => Cologne: 457066m
	Distance of Nuremberg from the 
	shortest path on earth
	(geodesic) between Munich and Cologne: 78745
	Point on Munich=>Cologne: 48.9467 10.3272 

