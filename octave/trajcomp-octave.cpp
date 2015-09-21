#include <octave/oct.h>
#include "../trajcomp.hpp"

#include <iostream>
#include <fstream>
using namespace std;

#define UNUSED(x) (void)(x)


// should put     OCTAVE_QUIT; into working loops to react on STRG-C
// http://fossies.org/windows/misc/octave-3.8.1.tar.gz/octave-3.8.1/doc/interpreter/external.texi#Input-Parameter-Checking-in-Oct_002dFiles

#define GLOBAL_BOUND_CHECK(x) { if ( x < 0 || x >= (int) trajectories.size()) \
			{ \
				octave_stdout << "bound error" << std::endl;\
				return octave_value_list(); \
			} }

typedef trajcomp::trajectory<double> trajectory;
trajcomp::default_segment_distance<trajectory::value_type> dseg;

std::vector<trajectory> trajectories;

// INIT
DEFUN_DLD(trajcomp_init, args,nargout,
          "Initialize Trajcomp Library")
{
	UNUSED(args);
	UNUSED(nargout);
   mlock(); // just to keep mlock consistently on.
            // I expect this method to be used before calling init_bgl
   trajectories.clear();
   return octave_value_list();
}


// load

DEFUN_DLD(trajcomp_load, args, nargout,
          "Load a trajectory")
{
	UNUSED(nargout);
    // Method not properly called.
    if (args.length() != 1) {
        octave_stdout << "Usage: trajcomp_load (A)\n";
        return octave_value_list();
    }

	std::string filename = args(0).string_value();
    trajectory t;
    t.load(filename);
			  	
	// Store the state
    trajectories.push_back(t);
    int index = trajectories.size()-1;
    return octave_value(index);
}


DEFUN_DLD(trajcomp_save, args, nargout,
          "Save a trajectory")
{
 	UNUSED(nargout);
    // Method not properly called.
    if (args.length() != 2 ) {
        octave_stdout << "Usage: @TODO\n";
        return octave_value_list();
    }

	int index = args(0).int_value();
	std::string filename = args(1).string_value();
    GLOBAL_BOUND_CHECK(index);
    
    trajectories[index].save(filename);

    return octave_value_list();
}


DEFUN_DLD(trajcomp_get, args, nargout,
          "Get a trajectory")
{
 	UNUSED(nargout);
    // Method not properly called.
    if (args.length() != 1) {
        octave_stdout << "Usage: @TODO\n";
        return octave_value_list();
    }

	int index = args(0).int_value();
    GLOBAL_BOUND_CHECK(index);
 
	Matrix A(trajectories[index].size(), trajectories[index].dimension);
    for (size_t i = 0; i < trajectories[index].size(); i++)
      for (size_t j=0; j <  trajectories[index].dimension; j++)
        A(i,j) = trajectories[index][i][j];

    return octave_value(A);
}


DEFUN_DLD(trajcomp_douglas_peucker, args, nargout,
          "Return a DP trajectory")
{
 	UNUSED(nargout);
    // Method not properly called.
    if (args.length() != 2) {
        octave_stdout << "Usage: @TODO\n";
        return octave_value_list();
    }

	int index = args(0).int_value();
	double epsilon = args(1).double_value();
    GLOBAL_BOUND_CHECK(index);
 
	trajectory t = trajcomp::douglas_peucker(trajectories[index],epsilon,dseg);
    trajectories.push_back(t);  
	unsigned int i = trajectories.size()-1;
    return octave_value(i);
}
