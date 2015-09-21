wgs84\_element\_distance: The distance between WGS84 coordinates
===================
It is not easy to correctly calculate the distance between WGS84 coordinates due to complexitites with the geometry of the WGS 84 ellipsoid. Our implementation of the distance (as well as the
point segment distance) are based on the C++ variant of the beautiful [GeographicLib library](http://geographiclib.sourceforge.net/) available at [http://geographiclib.sourceforge.net/](http://geographiclib.sourceforge.net/)

The problem of calculating the distance is tightly bound to the [Inverse Geodetic Problem](http://en.wikipedia.org/wiki/Geodesy), which we actually solve with the help of the aforementioned library. 

Dependency
-------------
If you write code that includes (and instatiates) the header trajcomp_geo.hpp, which defines this and some other geodetic functions, you have to install GeographicLib in the default place, such that the headers can be located. For Debian (and similar) Linux distributions, this is as simple as

	sudo apt-get install libgeographiclib-dev

putting the headers into the right place. Futhermore, GeographicLib is not a header-only library and you have to link with it using `-lGeographic` for g++


Quick Reference
-------------------
	// Class template
	template<class sample_type> 
	class wgs84_element_distance: public element_distance<sample_type,double>
	
	// Typical Use
	trajcomp::wgs84_element_distance<std::vector<double>> d;


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

