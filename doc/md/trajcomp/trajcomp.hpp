#ifndef TRAJCOMP_HPP_INC
#define TRAJCOMP_HPP_INC

#include<vector>
#include<iostream>
#include<stdexcept>
#include<string>
#include<sstream>
#include<fstream>
#include<limits>
#include<sys/time.h>
#include<algorithm>
#include<stdint.h>
#include<cmath>
#include<iomanip>

//#include<boost/filesystem.hpp>

#ifdef TRAJCOMP_SERIALIZATION
	// forward declaration for friendship
	namespace boost {
	namespace serialization {
		class access;
	}
	}
#endif


#include<math.h> // for sqrt
#define COL0 "\033[G"
#define CTLK "\033[K"
using namespace std;


//#define abs(x) (x<0?-x:x)

namespace trajcomp
{
	#define MIN(x,y) ((x<y)?x:y)
	#define MAX(x,y) ((x>y)?x:y)
	namespace tools{
		template<class datatype, char delim=' '>
		std::string make_string2(datatype data)
		{
			std::stringstream ss;
			typename datatype::iterator it;
			ss << std::endl;
			for(it=data.begin(); it != data.end(); it++)
			{
			  typename datatype::value_type::iterator it2;
			  for (it2 = (*it).begin();
				   it2 != (*it).end(); it2++){

					ss << "\t" << *it2 << delim;
				}ss  << std::endl;
			}
				return ss.str();
		}

		template<class datatype, char delim=' '>
		std::string make_string(datatype data)
		{
			std::stringstream ss;
			typename datatype::iterator it;
			for(it=data.begin(); it != data.end(); it++)
				ss <<  *it << delim;
			return ss.str();
		}

		void progress(long pos, long size,std::string what="")
		{
			std::cout << COL0 << "Progress("<<what<<"): " << pos << "/" << (double)size << "(" << 100*(double)pos/(double)size << "%)" << CTLK;
		}
		uint64_t ticks(void)
		{
			struct timeval tv;
			gettimeofday(&tv, 0);
			return uint64_t( tv.tv_sec ) * 1000 + tv.tv_usec / 1000;
		}
		class tictoc
		{
			public:
			uint64_t t;
			void tic() {t = ticks();};
			uint64_t toc(){return ticks() - t; };
		};

		template<class datatype>
		void matrix_resize(std::vector<std::vector< datatype> > &m, int r, int c)
		{
			m.resize(r);
			for(size_t i=0; i <r; i++)
			  m[i].resize(c);
		}
		template<class datatype>
		void matrix_resize(std::vector<std::vector<std::vector< datatype > > > &m, int r, int c, int d)
		{
			m.resize(r);
			for(size_t i=0; i <r; i++)
			{
			  m[i].resize(c);
			  for (size_t j=0; j < c; j++)
			    m[i][j].resize(d);
		    }

		}

		//------------------- BEGIN JANS IMPL ----------------------------
		double toRadiant(double val) {
			return val * (M_PI/180);
		}

		double toDegree(double val) {
			return val * (180/M_PI);
		}

		//-----------------START DENSITY SEGMENTATION----------------
		// Haversine distance between two coordinates (in km)360
		double haversine1(double lat1, double lat2, double lon1, double lon2) {
			const double earthRadius = 6378137;

			//cout << "Lat1: " << lat1 << endl;
			//cout << "Lon1: " << lon1 << endl;
			//cout << endl;
			double radLat1 = toRadiant(lat1);
			double radLat2 = toRadiant(lat2);

			//cout << "Lat1: " << radLat1 << endl;
			//cout << "Lon1: " << radiant(lon1) << endl;
			double dlon = toRadiant((lon2 - lon1));
			double dlat = toRadiant((lat2 - lat1));

			double a = pow(sin(dlat/2.0), 2) + cos(radLat1)*cos(radLat2)*pow(sin(dlon/2.0), 2);
			double c = 2 * atan2(sqrt(a), sqrt(1-a));
			//cout << "distance: " << earthRadius * c << endl;
			return earthRadius * c;
		}

		double calcBearing(double lat1, double lat2, double lon1, double lon2) {
			lat1 = toRadiant(lat1);
			lat2 = toRadiant(lat2);
			double dlon = toRadiant((lon2 - lon1));

			double X = sin(dlon) * cos(lat2);
			double Y = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon);
			double bearing = fmod((atan2(X,Y) * (180/M_PI) + 360.0), 360);
			cout << "bearing " << bearing << endl;
			return bearing;
		}

		double calc_orientation(double lon1, double lon2, double lat1, double lat2)
		{
			lon1 = trajcomp::tools::toRadiant(lon1);
			lon2 = trajcomp::tools::toRadiant(lon2);
			lat1 = trajcomp::tools::toRadiant(lat1);
			lat2 = trajcomp::tools::toRadiant(lat2);
			//http://www.sunearthtools.com/tools/distance.php#txtDist_3
			double delta_lon = lon2-lon1;

			if(abs(delta_lon) > M_PI) {
				if(delta_lon > 0.0) {
					delta_lon = -(2.0 * M_PI - delta_lon);
				} else {
					delta_lon = (2.0 * M_PI -delta_lon);
				}
			}

			double delta_phi = log(tan((lat2/2) + (M_PI/4)) / tan((lat1/2) + (M_PI/4)));

			double orientation = atan2(delta_lon , delta_phi);
			//orientation = fmod(trajcomp::tools::toDegree(orientation)+360.0, 360.0);
			return trajcomp::tools::toDegree(orientation);
		}

		/*
		* Berechnet den Geschwindigkeitsvektor zwischen zwei Punkten
		* returns vector<velocity, orientation>
		*/
		vector<double> calc_velocity_vector(double lon1, double lon2, double lat1, double lat2, double t1, double t2)
		{


			vector<double> velocity_vector;
			//velocity
			//cout << "P1: " << lon1 << "|" << lat1 << "|" << t1 << " - P2: " << lon2 << "|" << lat2 << "|" << t2 << endl;
			double mps = trajcomp::tools::haversine1(lat1, lat2, lon1, lon2) / ((t2 - t1) / 1000.0);
			//cout << "velocity: " << mps << endl;

			velocity_vector.push_back(mps);

			double orientation = trajcomp::tools::calc_orientation(lon1, lon2, lat1, lat2);

			velocity_vector.push_back(orientation);
			return velocity_vector;
		}

		template<class ttime>
		typename ttime::value_type calc_min_time(ttime &times)
		{
			typename ttime::value_type min_time = times[1] - times[0];
			for(size_t i = 1; i < times.size()-1; i++) {
				typename ttime::value_type cur_time = times[i+1] - times[i];
				if(cur_time < min_time) {
					min_time = cur_time;
				}
			}
			std::cout << "MIN TIME: " << min_time/1000 << endl;
			return min_time/1000;
		}

		/*
		*  Get Noise points by comparing all values of original and segmented
		*  trajectories. They have same ordering, so when a value of segmented
		*  is not in original, it is noise.
		*/
		template<class ttraj>
		ttraj getNoisePoints(ttraj &original, ttraj &segmented)
		{
			float eps = 0.00000001;
			ttraj noise(2);
			for(size_t i = 0; i < segmented.size(); i++) {
				if((fabs(segmented[i][0] - original[i][0]) < eps) &&
					(fabs(segmented[i][1] - original[i][1]) < eps)){
						continue;
				} else {
					noise.push_back({original[i][0], original[i][1]});
					original.erase(original.begin()+i);
				}
			}
			return noise;
		}


		/**
		* Read data from file, resample, write to given location and return trajectory
		* Source of algorithms: http://www.movable-type.co.uk/scripts/latlong.html
		* @param original Original trajectory
		* @param res Resulting trajectory after resampling
		* @param location_name The location to which resulting trajectory will be written
		* @param h Given distance to resample after
		*/
		template<class ttraj, class ttime>
		size_t resample_trajectory(ttraj &original, ttime &o_time, ttraj &res,
			ttime &r_time, std::string location_name, std::string dataName, double check_meters)
		{
			const double earthRadius = 6378137;
			typename ttime::value_type min_time = calc_min_time(o_time);
			//Calculate velocity vector between two points
			//Initialise result writers

			std::stringstream dataLocation;
			std::stringstream timeLocation;
			dataLocation << location_name << dataName << ".csv";
			timeLocation << location_name << dataName << "Time.csv";
			cout << dataLocation.str() << "   " << timeLocation.str() << endl;
			std::ofstream traj_out(dataLocation.str().c_str());
			traj_out.precision(10);
			std::ofstream time_out(timeLocation.str().c_str());
			time_out.precision(13);
			if (!traj_out.is_open() || !time_out.is_open())
			{
				throw(std::runtime_error("Could not open file."));
			}
			res.clear();
			r_time.clear();
			typename ttraj::value_type p1;
			typename ttraj::value_type p2;
			typename ttime::value_type pt1;
			typename ttime::value_type pt2;

			res.push_back(original[0]);
			r_time.push_back(o_time[0]);
			traj_out << original[0][0] << " " << original[0][1] << std::endl;
			time_out << o_time[0]<< std::endl;

			for(size_t i = 0; i < original.size()-1; i++) {
				p1.clear();
				p1 = original[i];
				pt1 = o_time[i];

				p2.clear();
				pt2 = o_time[i+1];
				p2 = original[i+1];

				double d = haversine1(p1[0], p2[0], p1[1], p2[1]);
				double velocity = d / ((pt2 - pt1 )/1000);
				double cur_progress = check_meters;
				double time_diff = 0.0;
				//Advance check_meters in every step until reaching distance to next point
				while(cur_progress < d) {
					//Check if time difference is greater than the minimal time difference
					if((time_diff = cur_progress / velocity) > min_time && ((d - cur_progress) / velocity) > min_time) {
						cout << "Time diff: " << time_diff << endl;
						//calculate new point
						double bearing = toRadiant(calcBearing(p1[0], p2[0], p1[1], p2[1]));

						double angular_distance = cur_progress/earthRadius;
						double radLat = toRadiant(p1[0]);
						double radLon = toRadiant(p1[1]);

						double lat_new = asin(sin(radLat) * cos(angular_distance) +
							cos(radLat) * sin(angular_distance) * cos(bearing));

						double lon_new = radLon + atan2((sin(bearing)*sin(angular_distance)*cos(radLat)),
							(cos(angular_distance)-sin(radLat) * sin(lat_new)));
						lat_new = toDegree(lat_new);
						lon_new = toDegree(lon_new);
						p1[0] = lat_new;
						p1[1] = lon_new;
						//Seconds to Milliseconds
						pt1 = pt1 + (time_diff*1000);
						//Push new point to results
						res.push_back(p1);
						r_time.push_back(pt1);
						traj_out << p1[0] << " " << p1[1] << std::endl;
						time_out << pt1<< std::endl;

						d = haversine1(p1[0], p2[0], p1[1], p2[1]);
						cur_progress = 0;
					}
					cur_progress += check_meters;
				}
				//Push last point to results
				res.push_back(p2);
				r_time.push_back(pt2);
				traj_out << p2[0] << " " << p2[1] << std::endl;
				time_out << pt2<< std::endl;
			}
			traj_out.close();
			time_out.close();
			std::cout << "Done resampling " << endl;
			return res.size();
		}

		//-------------------- END JANS IMPL--------------------------------

	} // tools

template <class ElementType, class DistanceType=double>
class element_distance
{
public:
   double operator() (void ){return 0;};
   virtual double operator()(ElementType u, ElementType v )=0;
};



template <class ElementType, class DistanceType=double>
class element_segment_distance
{
public:
   virtual DistanceType operator()(ElementType s1, ElementType s2, ElementType p)=0;
};

/*Default distances for containers of equal length*/

template<class sample_type>
class default_element_distance: public element_distance<sample_type,double>
{
	public:
	double operator() (void ){return 0;};
   virtual double operator()(sample_type u, sample_type v )
  {
	  if (u.size() != v.size())
		return std::numeric_limits<double>::quiet_NaN();
      double sum=0;
      for (size_t i=0; i < u.size(); i++)
        sum += (u[i]-v[i])*(u[i]-v[i]);

	return sqrt(sum);
  }
};

template<class sample_type>
class default_element_distance_squared: public element_distance<sample_type,double>
{
	public:
   virtual double operator()(sample_type u, sample_type v )
  {
	  if (u.size() != v.size())
		return std::numeric_limits<double>::quiet_NaN();
      double sum=0;
      for (size_t i=0; i < u.size(); i++)
        sum += (u[i]-v[i])*(u[i]-v[i]);

	return sum;
  }
};


template<class sample_type>
class default_segment_distance: public element_segment_distance<sample_type,double>
{
	public:
   virtual double operator()(sample_type u,sample_type v, sample_type p)
  {
	  // dimensions equal?
	  if (u.size() != v.size() || v.size() != p.size())
		return std::numeric_limits<double>::quiet_NaN();

	  // the length of the segment
	  default_element_distance_squared<sample_type> d2;

	  double l = d2(u,v);
#ifdef DEBUG_default_segment_distance
	  std::cout << "l=" << l << endl;
#endif

	  if (fabs(l) < 1E-12)
			return sqrt(d2(p,v));

		double t = 0;
		for (size_t i=0; i <  u.size(); i++)
		   //t += ((u[i]-v[i])*(p[i]-v[i])); @TODO:Remove?
		   t += ((v[i]-u[i])*(p[i]-u[i]));
		t/=l;
#ifdef DEBUG_default_segment_distance
	std::cout << "t= " << t << endl;
#endif

		if (t<0) return sqrt(d2(p,u));
		if (t>1) return sqrt(d2(p,v));

		sample_type proj;
		for (size_t i=0; i <  u.size(); i++)
		   //proj.push_back(v[i]+t*(u[i]-v[i]));
		   proj.push_back(u[i]+t*(v[i]-u[i]));
		double ret = sqrt(d2(p,proj));
#ifdef DEBUG_default_segment_distance
		std::cout << "(" <<
				tools::make_string(u) << ";" <<
				tools::make_string(v) << ";" <<
				tools::make_string(p) << ") ==>";
		std::cout << tools::make_string(proj) << "=" << ret <<  std::endl;
#endif
		return ret;
  }
};





	// A trajectory is extended from a vector of elements of a dimension
	template <typename ValueType>
	class trajectory : public std::vector<std::vector<ValueType>>
	{

		public:
		unsigned int dimension;
		typedef std::vector<std::vector<ValueType>> Base;
		typedef std::vector<ValueType> ElementType;

		trajectory():dimension(0)
		{};
		trajectory(size_t dim):dimension(dim)
		{};
		void summary()
		{
			std::cout << "Dimension = " << dimension << std::endl;
			std::cout << "Elements  = " << this->size() << std::endl;
		}

		void save(std::string filename)
		{
			std::ofstream of(filename);
			if (of)
			{
				for (typename Base::iterator it=this->begin(); it != this->end(); it++)
				{
					std::string s = tools::make_string(*it);
					of << s << std::endl;

				}
			}else
			{
				throw(std::runtime_error("Error: Can't save to file."));
			}

		}


		void load(std::string filename)
		{
			this->clear();
			std::ifstream infile(filename);
			if (!infile)
			{
				throw(std::runtime_error("File not found."));
			}
			std::string line;
			while (std::getline(infile, line))
			{
				ElementType l;
				std::stringstream iss(line);
				ValueType v;
				while (iss >> v)
					l.push_back(v);
				this->dimension=l.size();
				this->push_back(l);
			}
		}


		void dump()
		{
			for (typename Base::iterator it=this->begin(); it != this->end(); it++)
		   {
			  std::string s = tools::make_string(*it);
			  std::cout << s << std::endl;
		   }
		}

		void push_back (const typename Base::value_type& val)
		{
			/*if(val.size() != dimension)
			 throw(std::runtime_error("Error: Can't add element of wrong dimension to trajectory."));*/
			Base::push_back(val);

		}
		void push_back (typename Base::value_type && val)
		{
			/*if(val.size() != dimension)
			 throw(std::runtime_error("Error: Can't add element of wrong dimension to trajectory."));*/
			Base::push_back(val);
		}

		/*SERIALIZATION*/

#ifdef TRAJCOMP_SERIALIZATION
		 // Allow serialization to access non-public data members.
		friend class boost::serialization::access;

		template<typename Archive>
		void serialize(Archive& ar, const unsigned version) {
				ar & d1_ & d2_;  // Simply serialize the data members of Obj
		}
#endif

	};


//3.2.1. Uniform Select.
// 	Extracts each k-th point

template<class TrajectoryType>
TrajectoryType uniform_select(TrajectoryType &t,size_t Step=15,bool withEndPoint=true)
{
	TrajectoryType ret(t.dimension);
	typename TrajectoryType::iterator it;
	size_t k=0;
    for (it = t.begin(); it != t.end(); it++)
    {
		if (k % Step == 0)
			ret.push_back(*it);
		k++;
	}

	if(withEndPoint)
	  if (((t.size()-1) % Step) != 0)
	  {
	    ret.push_back(t.back());
	}
	return ret;
}

/*
template<class TrajectoryType, typename functional>
TrajectoryType middle_distance_interpolation(TrajectoryType &t,double d_max,
											functional f)
{
	TrajectoryType ret(t.dimension);
	typename TrajectoryType::iterator it;

	first = t.begin();
	auto second = first;
	second ++;
    for (; second != t.end();)
    {
		double seglen = f(*first, *second);

		first ++;
		second ++;
	}

	return ret;
}
*/

//---------------------- START JANS IMPL-----------------------------------
template <class TrajectoryType, class DistanceType=double, class PtsType=int>
class DBSCAN_segmentation_impl
{
	public:
		std::string exec(const char* cmd)
		{
			FILE* pipe = popen(cmd, "r");
		    if (!pipe) return "ERROR";
		    char buffer[128];
		    std::string result = "";
		    while(!feof(pipe)) {
		    	if(fgets(buffer, 128, pipe) != NULL)
		    		result += buffer;
		    }
		    pclose(pipe);
		    return result;
		}
		double get_mean_in_cluster(vector<vector<double>> &coords)
		{
			double mean = 0.0;
			for(size_t i = 0; i < coords.size(); i++) {
				mean += coords[i][1];
			}
			mean = mean/coords.size();
			size_t meanID = coords.size()-1;
			for(size_t j = 0; j < coords.size()-1; j++) {
				//std::cout << meanID << " ID; Current diff to mean: " << fabs(coords[meanID][1] - mean);
				//std::cout << " " << coords[j][0] << " ID; Next diff to mean: " << fabs(coords[j][1] - mean) << endl;
				if(fabs(coords[j][1] - mean) < fabs(coords[meanID][1] - mean)) {
				//	std::cout << "FOUND SMALLER " << coords[j][0] << endl << endl;
					meanID =  j;
				}
			}
			//std::cout << "Found meanest: " << coords[meanID][0] << endl;
			return coords[meanID][0];
		}

		bool sortByID(const vector<double> &a, const vector<double> &b){return a[0] < b[0];}

		TrajectoryType operator()(TrajectoryType &traj, DistanceType epsilon,
			PtsType minPts, const std::string &dataName, const std::string &pathToData,
			const std::string path_to_elki)
		{
			//Schritt 1: sample data (Evtl auch davor in der main aufrufen)
			std::cout << "start dbscan..." << std::endl;
			//Schritt 2: Rufe DBSCAN auf übergebenen Datenpfad auf
			std::stringstream cmd;
			cout <<  "data name: " << dataName << endl;
			cout << "path to data: " << pathToData << endl;
			cmd << "java -jar " << path_to_elki << " KDDCLIApplication -dbc.in " << pathToData <<
			" -algorithm clustering.DBSCAN -dbscan.epsilon " << epsilon << " -dbscan.minpts " <<
			minPts << " -resulthandler ResultWriter -out clustering/" << dataName;
			const std::string& tmp = cmd.str();
			const char* p = tmp.c_str();
			string result = exec(p);

			//cout << result << endl;

			//TODO Überprüfen, ob bereits in trajektorie vorhanden
			//Erster Punkt von Trajektorie muss im Ergebnis sein
			vector<size_t> resIDs;
			resIDs.push_back(0);
			TrajectoryType res(2);
			res.clear();

			//{ID, lat+lon}
			vector<vector<double> > coords(2);

			//Gehe gefundene Cluster in Ordner durch
			int i = 0;
			while(true)
			{
				std::stringstream path_name;
				path_name << "clustering/" << dataName << "/cluster_" << i << ".txt";
				std::ifstream file(path_name.str().c_str());
				if(!file) {
					cout << "No more clusters found... " << endl;
					break;
				}
			//cout << "Found file " << endl;
				coords.clear();
				string line;
				while(getline(file, line)) {
					size_t found = line.find("#");

					if(found != string::npos){
						continue;
					}
					std::istringstream iss(line);
					string IDString;
					double lat, lon;
					if(iss >> IDString >> lat >> lon) {
					} else {
						throw(std::runtime_error("Error: Can't convert values in cluster."));
					}
					IDString = IDString.substr(IDString.find("=")+1, string::npos);
					iss.str("");
 					iss.clear();
					iss.str(IDString);
					double ID;
					iss >> ID;
					//std::cout << "CLUSTERID: " << ID << std::endl;
					//std::cout << lat << " : " << lon << std::endl;
					coords.push_back({ID,lat+lon});
				}
				//Schritt 3: Suche repräsentierbare Punkte in Cluster
				size_t meanestID = get_mean_in_cluster(coords);
				resIDs.push_back(meanestID);
		    	file.close();
				i++;
			}

			//Letzter Punkt in Trajektorie muss im Ergebnis sein
			if(std::find(std::begin(resIDs), std::end(resIDs), traj.size()-1) != std::end(resIDs)) {
				resIDs.push_back(traj.size()-1);
			}
			//Delete cluster data
			std::sort(resIDs.begin(), resIDs.end());
			for(int i = 0; i < resIDs.size(); i++) {
				res.push_back(traj[resIDs[i]]);
			}
			return res;

		}

};
/*
* @param traj the trajectory to be segmented loaded as trajectory datatype
* @param epsilon user defined epsilon
* @param minPts user defined minimal points in cluster
* @param dataName name of data, clusters will be save in a folder with this name
* @param pathToData path to data that should be segmented
* @param path_to_elki path to elki.jar
*/
template <class TrajectoryType>
TrajectoryType DBSCAN_segmentation(TrajectoryType &traj, double epsilon,
	int minPts, string dataName, string pathToData, string path_to_elki)
{
	DBSCAN_segmentation_impl<TrajectoryType, double, int> dbs;
	//boost::filesystem::remove_all("clustering/");
	system("exec rm -r clustering/*");
	TrajectoryType ret = dbs(traj,epsilon, minPts, dataName, pathToData, path_to_elki);

	//TODO: What to do with noise points?

	return ret;
}

//--------------------START THRESHOLD SAMPLING---------------------------


template <class TrajectoryType, class TimesType, class ThresholdType=double>
class threshold_sampling_impl
{
	public:

		TrajectoryType operator()(TrajectoryType &traj, TimesType &times,
			ThresholdType velocity_thresh, ThresholdType orientation_thresh)
		{
			std::cout << std::setprecision(10);
			TrajectoryType result(2);
			result.clear();
			TimesType resultTimes(1);
			resultTimes.clear();
			result.push_back({traj[0][0], traj[0][1]});

			resultTimes.push_back(times[0]);
			result.push_back({traj[1][0], traj[1][1]});
			resultTimes.push_back(times[1]);

			//calculate safe area
			for(size_t i = 2; i<traj.size(); i++) {
				size_t resultSize = result.size();

				//START Calculate upper and lower angle and radius for trajectory and sample based last two points.

				//For velocity vector of trajectory, take last 2 elements in original trajectory
				vector<double> velocity_vector_traj = trajcomp::tools::calc_velocity_vector(traj[i-2][1], traj[i-1][1],
													traj[i-2][0], traj[i-1][0], times[i-2], times[i-1]);
				//For velocity vector of sample, take last 2 elements in result trajectory
				vector<double> velocity_vector_sample = trajcomp::tools::calc_velocity_vector(result[resultSize-2][1], result[resultSize-1][1], result[resultSize-2][0],
												result[resultSize-1][0], resultTimes[resultSize-2], resultTimes[resultSize-1]);

				//Calculate the radius for trajectory and sample based using the calculated velocity vectors
				double upper_radius_traj = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_traj[0]*(1 + velocity_thresh));
				double lower_radius_traj = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_traj[0]*(1 - velocity_thresh));

				double upper_radius_sample = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_sample[0]*(1 + velocity_thresh));
				double lower_radius_sample = ((times[i] - times[i-1]) / 1000.0) *
												(velocity_vector_sample[0]*(1 - velocity_thresh));

				cout << "Upper radius traj: " << upper_radius_traj << " lower: " << lower_radius_traj << endl;
				cout << "Uppe radius sample: " << upper_radius_sample << " lower: " << lower_radius_sample << endl;
				//Calculate orientation angle for both from given slopes
				double upper_slope_traj = velocity_vector_traj[1] + orientation_thresh;
				double lower_slope_traj = velocity_vector_traj[1] - orientation_thresh;


				double upper_slope_sample = velocity_vector_sample[1] + orientation_thresh;
				double lower_slope_sample = velocity_vector_sample[1] - orientation_thresh;
				cout << "Upper slope traj: " << upper_slope_traj << " lower: " << lower_slope_traj << endl;
				cout << "Uppe slope sample: " << upper_slope_sample << " lower: " << lower_slope_sample << endl;
				//START Calculate angle and radius for point to insert
				vector<double> velocity_vector_insert = trajcomp::tools::calc_velocity_vector(traj[i-1][1], traj[i][1],
													traj[i-1][0], traj[i][0], times[i-1], times[i]);

				double radius_insert = ((times[i] - times[i-1]) / 1000.0) *
												velocity_vector_insert[0];
				cout << "Radius insert: " << radius_insert <<endl;
				double angle_insert = velocity_vector_insert[1];
				cout << "angle insert: " << angle_insert << endl;

				//Check if new points angle and radius is inside of the calculated
				bool isInRadius = false;
				bool isInAngle = false;
				if((radius_insert >= lower_radius_traj) && (radius_insert >= lower_radius_sample)
						&& (radius_insert <= upper_radius_traj) && radius_insert <= upper_radius_sample) {
					isInRadius = true;
				}

				if((angle_insert <= upper_slope_traj) && (angle_insert <= upper_slope_sample)
						&& (angle_insert >= lower_slope_sample) && (angle_insert >= lower_slope_traj)) {
					isInAngle = true;
				}
				if(!isInAngle || !isInRadius) {
					result.push_back({traj[i][0], traj[i][1]});
					resultTimes.push_back(times[i]);
				}
			}

			cout << "Original size: " << traj.size() << " threshold sampling result size: " << result.size() << endl;

			return result;
		}
};

/*
*
*/
template<class TrajectoryType, class TimesType>
TrajectoryType threshold_sampling(TrajectoryType &traj, TimesType &times, double velocity_thresh, double orientation_thresh)
{
	threshold_sampling_impl<TrajectoryType, TimesType, double> ts;

	TrajectoryType ret = ts(traj, times, velocity_thresh, orientation_thresh);
	return ret;
}


//--------------------END JANS IMPL----------------------------------------

//3.2.2. Douglas Peucker Algorithm.
// Implemented using a complex class to encapsulate methods for
// recursion and a default function template

template <class TrajectoryType, class DistanceType=double>
class douglas_peucker_impl
{
	public:
	std::vector<unsigned char> booleanMap; // @REMARK: could also use bool specialization
	element_segment_distance<typename TrajectoryType::value_type, DistanceType> *d;
public:
   douglas_peucker_impl():d(NULL){};
DistanceType largest_contribution(TrajectoryType &u, size_t istart, size_t iend,size_t &the_index)
   {
	// find the point which contributes most

	DistanceType the_contribution = -1;
	typename TrajectoryType::value_type s = u[istart];
	typename TrajectoryType::value_type e = u[iend];

	for (size_t i=istart+1; i < iend; i++)
	{

		DistanceType contribution = (*d)(s,e,u[i]);
#ifdef DEBUG_largest_contribution
		std::cout << i << "~" << contribution << std::endl;
#endif
		if (contribution > the_contribution)
		{
			the_contribution = contribution;
			the_index = i;
		}
	}

	return the_contribution;
   }

   void simplify(TrajectoryType &u, size_t istart, size_t iend, DistanceType max_error)
   {
	   // simplify segment and recurse
#ifdef DEBUG_douglas_peucker_simplify
	   std::cout << istart << "," << iend << std::endl;
#endif


	   if ((istart -iend) < 2) return; // there is nothing I can add.


	   size_t the_index=0;
	   DistanceType the_contribution = largest_contribution(u,istart,iend,the_index);
	   if (the_contribution > max_error)
	   {
#ifdef DEBUG_douglas_peucker_simplify
		   std::cout << "DP: Insert " << the_index << std::endl;
#endif
		   booleanMap[the_index] = 1;
		   simplify(u,istart,the_index,max_error);
		   simplify(u,the_index,iend,max_error);
	   }
   }

   TrajectoryType operator()(TrajectoryType u, DistanceType max_error,
            element_segment_distance<typename TrajectoryType::value_type,DistanceType> &d)
  {
    this->d = &d;

	booleanMap.resize(u.size());
	for (size_t i=0; i < u.size(); i++)
		booleanMap[i] = 0;
	booleanMap[0] = booleanMap[u.size()-1] = 1;

	simplify(u,0,u.size()-1,max_error);
	// now build the simplified trajectory
	TrajectoryType ret;
	for (size_t i=0; i < booleanMap.size(); i++)
	  if (booleanMap[i])
	     ret.push_back(u[i]);

    return ret;
  }
};



template<class TrajectoryType>
TrajectoryType douglas_peucker(TrajectoryType &t, double epsilon,
element_segment_distance<typename TrajectoryType::value_type,double> &d)
{
	douglas_peucker_impl<TrajectoryType, double> dp;
	TrajectoryType ret = dp(t,epsilon,d);
	return ret;
}


template<class TrajectoryType>
TrajectoryType douglas_peucker(TrajectoryType &t, double epsilon)
{
	douglas_peucker_impl<TrajectoryType, double> dp;
	default_segment_distance<typename TrajectoryType::value_type> d;
	TrajectoryType ret = dp(t,epsilon,d);
	return ret;
}



//3.2.3. Bellmans Algorithm.
//3.2.4. Bottom up.
//3.2.5. Reservoir Sampling.
//3.2.6. Sliding Window Algorithm.

//////////////////// ELEMENTARY DISTANCE ALGORITHMS

//// POINT TO TRAJECTORY DISTANCE

template<class TrajectoryType>
typename TrajectoryType::ElementType::value_type point_to_trajectory_point(typename TrajectoryType::ElementType &p, TrajectoryType &t)
{
	default_element_distance<typename TrajectoryType::ElementType> d;
	typename TrajectoryType::ElementType::value_type m = 0;
	bool first = true;

	for (typename TrajectoryType::iterator it=t.begin(); it != t.end(); it++)
	{
		typename TrajectoryType::ElementType::value_type v=d(p,*it);
		if (first || v<m)
		{
			   m=v;
			   first=false;
		}

	}
	return m;
}

template<class TrajectoryType>
typename TrajectoryType::ElementType::value_type point_to_trajectory(typename TrajectoryType::ElementType &p, TrajectoryType &t)
{
	default_segment_distance<typename TrajectoryType::ElementType> d;
	typename TrajectoryType::ElementType::value_type m = 0;

	bool first = true;
	for (size_t i=0; i < t.size()-1; i++)
	{
		typename TrajectoryType::ElementType::value_type v=d(t[i],t[i+1],p);
		if (first || v<m)
		{
			   m=v;
			   first=false;
		}

	}

	return m;
}




template<class TrajectoryType>
typename TrajectoryType::ElementType::value_type closest_pair(TrajectoryType &u, TrajectoryType &v)
{
	default_element_distance<typename TrajectoryType::ElementType> d;
	typename TrajectoryType::ElementType::value_type m = 0;
	bool first = true;

	for (typename TrajectoryType::iterator it1=u.begin(); it1 != u.end(); it1++)
	for (typename TrajectoryType::iterator it2=v.begin(); it2 != v.end(); it2++)
	{
		typename TrajectoryType::ElementType::value_type v=d(*it1,*it2);
		if (first || v<m)
		{
			   m=v;
			   first=false;
		}

	}
	return m;
}


template<class TrajectoryType>
class implicit_distance_matrix
{
	private:
	size_t cols;
	std::vector<typename TrajectoryType::ElementType::value_type> D;
	TrajectoryType &u;
	TrajectoryType &v;
	default_element_distance<typename TrajectoryType::ElementType> d;

	public:
		implicit_distance_matrix(TrajectoryType &theu,TrajectoryType &thev): u(theu), v(thev)
		{
			cols = v.size();
			D.resize(u.size()*v.size());
			for (size_t i=0; i < u.size() * v.size(); i++)
			  D[i] = -1;
		}

	typename TrajectoryType::ElementType::value_type  operator ()(size_t i, size_t j)
	{
		size_t index = i*cols + j;
		if (D[index] < 0)
		  D[index] = d(u[i],v[j]);
		return D[index];
	}
};




template<class TrajectoryType>
typename TrajectoryType::ElementType::value_type sum_of_pairs(TrajectoryType &u, TrajectoryType &v)
{
	if (u.size() != v.size())
	   throw(std::runtime_error("sum_of_pairs distance only defined for trajectories of equal size"));

	default_element_distance<typename TrajectoryType::ElementType> d;
	typename TrajectoryType::ElementType::value_type m = 0;
	typename TrajectoryType::iterator it1=u.begin();
	typename TrajectoryType::iterator it2=v.begin();

	while (it1 != u.end()) // equal length has been tested
	{
		typename TrajectoryType::ElementType::value_type v=d(*it1,*it2);
		m += v;
		it1++;
		it2++;

	}
	return m;

}






//5.1.1 Discrete Frechet Distance
template <class TrajectoryType, class DistanceType=double>
class discrete_frechet_distance_impl
{
	private:
	std::vector<std::vector<DistanceType > > CA;
	element_distance<typename TrajectoryType::ElementType,DistanceType> *d;
	protected:

	double c(int i,int j,TrajectoryType &t1, TrajectoryType &t2	)
	{
		if (CA[i][j]>-1)		// don't update
			return CA[i][j];
    if (i==0 && j==0){
        return (CA[i][j] = (*d)(t1[0],t2[0]));
    }
    if( i>0 && j==0){
        CA[i][j] = MAX(
			c(i-1,0,t1,t2),
			(*d)(t1[i],t2[0]));
        return CA[i][j];
    }
    if( i==0 && j>0){
        CA[i][j] = MAX(
			c(0,j-1,t1,t2),
			(*d)(t1[0],t2[j]));
        return CA[i][j];
    }
    if (i>0 && j>0)
    {
	 CA[i][j] = MAX(MIN(MIN(
			c(i-1,j  ,t1,t2),
		    c(i-1,j-1,t1,t2)),
			c(i  ,j-1,t1,t2)),
			(*d)(t1[i],t2[j]));
        return CA[i][j];
    }
    CA[i][j] = std::numeric_limits<double>::infinity();
    return CA[i][j];
}



public:
	discrete_frechet_distance_impl():d(NULL){};
    DistanceType operator()(TrajectoryType u, TrajectoryType v,
    element_distance<typename TrajectoryType::ElementType,DistanceType> &d)
  {
	  this->d = &d;
#ifdef DEBUG_DISCRETE_FRECHET
	 std::cout << "Comparing " << u.size() << " and " << v.size() << std::endl;
#endif
	DistanceType distance=0;
    if (u.size() == 0) return 0;
    for (size_t i=0; i < u.size(); i++)
    {
		CA.push_back(std::vector<DistanceType>());
       for (size_t j=0; j < v.size(); j++)
          CA[i].push_back(-1);
	  }
	  distance = c(u.size()-1,v.size()-1,u,v);
	return distance;
  }

};

template<class TrajectoryType>
double discrete_frechet(TrajectoryType &u, TrajectoryType &v)
{
	discrete_frechet_distance_impl<TrajectoryType, double> df;
	default_element_distance<typename TrajectoryType::ElementType> d;
	return df(u,v,d);
}



} // NAMESPACE


#endif
