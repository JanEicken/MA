#ifndef TRAJCOMP_GOOGLENOW_HPP
#define TRAJCOMP_GOOGLENOW_HPP

#include "trajcomp.hpp"
#include "trajcomp_files.hpp"
#include<fstream>

#include <sys/types.h>
#include <dirent.h>
#include <string.h>

#include<algorithm>


namespace trajcomp{


	// Dataset reader with time

/*
 * Line 1...6 are useless in this dataset, and can be ignored. Points are described in following lines, one for each line.
Field 1: Latitude in decimal degrees.
Field 2: Longitude in decimal degrees.
Field 3: All set to 0 for this dataset.
Field 4: Altitude in feet (-777 if not valid).
Field 5: Date - number of days (with fractional part) that have passed since 12/30/1899.
Field 6: Date as a string.
Field 7: Time as a string.
Note that field 5 and field 6&7 represent the same date/time in this dataset. You may use either of them.

 *
 * */




template<class TrajectoryType, size_t Lat=0,size_t  Lon=1>
class googleNow_reader
{
	public:
	std::string googleNow_base;
	size_t MAX_DIR;
	googleNow_reader(std::string base, size_t maxdir):googleNow_base(base),MAX_DIR(maxdir)
	{};
	googleNow_reader(std::string base):googleNow_base(base),MAX_DIR(200)
	{}

	TrajectoryType handle_file (const char *dir, const char *fname)
	{
		if (fname[0] == '.')
		{
			throw std::runtime_error("Invisible File");
		}
		char buffer[1024];

		strncpy(buffer,dir,1024);
		strncat(buffer,fname,1024);

		std::ifstream infile;
		infile.open(std::string(buffer));

		if (!infile)
		{
			throw(std::runtime_error("File not found."));
		}
		TrajectoryType traj;
		std::string line;
		while (std::getline(infile, line))
		{

			typename TrajectoryType::value_type element;
			element.resize(2);
			//std::cout << line << std::endl;

			std::stringstream iss(line);

			if (! (iss >> element[ Lat ]))
			  throw(std::runtime_error("Parsing Error at Lat"));
			if (! (iss >> element[Lon]))
			  throw(std::runtime_error("Parsing Error at Lon"));

			//std::cout << tools::make_string(element) << std::endl;

			 traj.push_back(element);
		}

    //traj.summary();
	return traj;
    }


	int load(std::vector<TrajectoryType> &db, bool isTime) //std::vector<TimeType> time
	{
		std::cout << "Loading from Base " << googleNow_base << std::endl;
		size_t i;
	    for (i=0; i < MAX_DIR ; i++)
       {
            char buffer[1024];
            if(isTime) {
                snprintf(buffer,1024,"%s/%03d/times/",googleNow_base.c_str(),i);
            } else {
                snprintf(buffer,1024,"%s/%03d/",googleNow_base.c_str(),i);
            }

#ifdef VERBOSE_GOOGLENOW_READER
            std::cout <<"Going to enumerate directory" << buffer << std::endl;
#endif

            /*DIR *dp;
			struct dirent *ep;*/

			struct dirent **namelist;
            int n;
            n = scandir(buffer, &namelist, NULL, alphasort);
            if (n < 0)
               perror("scandir");
           else {
			   for (size_t k=0; k < n; k++)
			   {
				   //intf("%s\n", namelist[k]->d_name);
				   TrajectoryType t;
				   try{
				   t = handle_file(buffer,namelist[k]->d_name);
			      }catch(std::runtime_error &e)
			      {
					  #ifdef VERBOSE_GOOGLENOW_READER
					  cout << "Caught " << e.what() << endl;
					  cout << "Skipping " << namelist[k]->d_name << endl;
					  #endif
					  continue;
				  }
					if (t.size() != 0)
					{
						db.push_back(t);
					}else{
						std::cout << "Warning: " << buffer << "led to empty trajectory." << std::endl;
					}
				   free(namelist[k]);
			   }


               free(namelist);
			}
		}
		return i;
	}


};

}

#endif
