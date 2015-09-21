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
