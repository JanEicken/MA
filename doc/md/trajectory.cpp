#ifdef _MSC_VER
#include <boost/config/compiler/visualc.hpp>
#endif
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <cassert>
#include <exception>
#include <time.h>

#include "trajcomp/trajcomp.hpp"
#include "trajcomp/trajcomp_files.hpp"

#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<algorithm>
#include<stdio.h>
#include<string>
#include<sstream>

using namespace std;

// specialize our types
//typedef trajcomp::trajectory<double> trajectory;
//START types used for all Segmentations
typedef vector<vector<double>> trajectory;
//END

//START Types used for Persistance based segmentation
//Typ curvatures: {curvature, isMaxima}
typedef vector< vector<int> > curvatures;
typedef vector<size_t> component;
//END

const double PI = 3.14159265;

double radiant(double val) {

	return val * (PI/180);
}

//-----------------START DENSITY SEGMENTATION----------------
// Haversine distance between two coordinates
double haversine(double lat1, double lat2, double lon1, double lon2) {
	const double earthRadius = 6371000;

	//cout << "Lat1: " << lat1 << endl;
	//cout << "Lon1: " << lon1 << endl;
	//cout << endl;
	double radLat1 = radiant(lat1);
	double radLat2 = radiant(lat2);

	//cout << "Lat1: " << radLat1 << endl;
	//cout << "Lon1: " << radiant(lon1) << endl;
	double dlon = radiant((lon2 - lon1));
	double dlat = radiant((lat2 - lat1));

	double a = pow(sin(dlat/2.0), 2) + cos(radLat1)*cos(radLat2)*pow(sin(dlon/2.0), 2);
	double c = 2 * atan2(sqrt(a), sqrt(1-a));
	//cout << "distance: " << earthRadius * c << endl;
	return earthRadius * c;
}
/*
* Idee hierbei:
* Entweder die Zeit hernehmen, um stay Orte zu finden, oder dies anhand der Häufung von Punkten tun
* Häufung von Punkten:
*	- Abstand von einem Punkt zum nächsten unterschreitet einen bestimmten Threshold. Wenn das so ist,
* 	betrachte Abstand zum wieder nächsten bis einer gefunden wird der drüber ist. Schreibe diesen als nächsten
*	Punkt in result trajectory
*	- Richtungswechsel sind egal!
*/
// Könnte auch so geschrieben werden, dass die Distanzfunktion mitübergeben wird.
template<class trajectory>
trajectory density_segmentation(trajectory &traj, double max_distance)
{

	trajectory res;
	res.push_back(traj[0]);
	for (size_t i =0; i < traj.size(); i++)
	{
		for(size_t j = i+1; j < traj.size(); j++)
		{
			double distance = haversine(traj[i][1], traj[j][1], traj[i][0], traj[j][0]);
			//std::cout << distance << std::endl;
			if(distance > max_distance) {
				res.push_back(traj[j]);
				i = j;
				break;
			} else {
				//entferne Knoten aus time_property
			}

		}
		//auto time = t[i];

	}
	res.push_back(traj[traj.size()-1]);
	return res;
}
//----------------------END DENSITY SEGMENTATION--------------------


//----------------------START ONLINE ALGORITHM------------------------------
/*
*Diese Funktion berechnet den Krümmungsradius an jedem Punkt.
* Der Abstand h sollte möglich klein gewählt werden
*
*/
template<class TrajectoryType, class DistanceType=double>
TrajectoryType calc_curvatur(TrajectoryType &traj, DistanceType h)
{
	TrajectoryType res = TrajectoryType(2);
	for(size_t i = 0; i<traj.size(); i++) {
		if(i == 0) {
			//Berechne nächsten Punkt mit Abstand h
			double d_next = sqrt(pow(traj[i+1][0]-traj[i][0],2) + pow(traj[i+1][1]-traj[i][1],2));
			vector<double> p_next(2);
			p_next.clear();

			if(traj[i+1][0] < traj[i][0]){
				p_next.push_back(traj[i][0] + traj[i][0]*(h/d_next));
			} else {
				p_next.push_back(traj[i][0] - traj[i][0]*(h/d_next));
			}
			if(traj[i+1][1] < traj[i][1]){
				p_next.push_back(traj[i][1] + traj[i][1]*(h/d_next));
			} else {
				p_next.push_back(traj[i][1] - traj[i][1]*(h/d_next));
			}
			//Default case
			//res.push_back({0.0,0.0});


		} else if (i == traj.size()-1) {
			//Berechne vorherigen Punkt mit Abstand h
			double d_prev = sqrt(pow(traj[i][0]-traj[i-1][0],2) + pow(traj[i][1]-traj[i-1][1],2));
			vector<double> p_prev(2);
			p_prev.clear();

			if(traj[i][0] < traj[i-1][0]){
				p_prev.push_back(traj[i][0] + traj[i][0]*(h/d_prev));
			} else {
				p_prev.push_back(traj[i][0] - traj[i][0]*(h/d_prev));
			}
			if(traj[i][1] < traj[i-1][1]){
				p_prev.push_back(traj[i][1] + traj[i][1]*(h/d_prev));
			} else {
				p_prev.push_back(traj[i][1] - traj[i][1]*(h/d_prev));
			}
			//Default case
			//res.push_back({0.0,0.0});

		} else {
			//Berechne zwei Punkte mit Abstand h zu gegenwärtigem Punkt
			double d_prev = sqrt(pow(traj[i][0]-traj[i-1][0],2) + pow(traj[i][1]-traj[i-1][1],2));
			double d_next = sqrt(pow(traj[i+1][0]-traj[i][0],2) + pow(traj[i+1][1]-traj[i][1],2));
			vector<double> p_prev(2);
			vector<double> p_next(2);
			p_prev.clear();
			p_next.clear();

			if(traj[i][0] < traj[i-1][0]){
				p_prev.push_back(traj[i][0] + traj[i][0]*(h/d_prev));
			} else {
				p_prev.push_back(traj[i][0] - traj[i][0]*(h/d_prev));
			}
			if(traj[i][1] < traj[i-1][1]){
				p_prev.push_back(traj[i][1] + traj[i][1]*(h/d_prev));
			} else {
				p_prev.push_back(traj[i][1] - traj[i][1]*(h/d_prev));
			}

			if(traj[i+1][0] < traj[i][0]){
				p_next.push_back(traj[i][0] + traj[i][0]*(h/d_next));
			} else {
				p_next.push_back(traj[i][0] - traj[i][0]*(h/d_next));
			}
			if(traj[i+1][1] < traj[i][1]){
				p_next.push_back(traj[i][1] + traj[i][1]*(h/d_next));
			} else {
				p_next.push_back(traj[i][1] - traj[i][1]*(h/d_next));
			}
			//double first_dev = (1/pow(h,2)) * ()

		}

	}

}

/*
* Berechnet lokale Minima und Maxima der gegebenen Trajektorie.
* Die Eingabe muss aus Krümmungsradien bestehen, minima und maxima werden in übergebenen vectoren gescpeichert
* In global_minima und global_maxima wird der Index der jeweiligen globalen extrema gespeichert
*
*/
template<class ElementType, class IndexType=vector<size_t>>
void calc_extrema(ElementType &curv, IndexType &minima)
{
	bool isClimbing = true;
	minima.clear();
	//Prüfe erstes Element der Trajektorie auf Minimum oder Maximum
	if(curv[0][0] < curv[1][0]) {
		minima.push_back(0);
		curv[0][2] = 1;
	} else {

		//curv[0][1] = 1;
		isClimbing = false;
	}
	//Prüfe innere Elemente auf Minima und Maxima
	for(size_t i = 1; i < curv.size(); i++) {
		//Last element
		if(i == curv.size()-1) {
			//cout << "last element: " << curv[i][0] << " climbing: " << isClimbing <<endl;
			if(isClimbing) {
				if(curv[i][0] > curv[i-1][0]) {
					curv[i][1] = 1;
				} else {
					minima.push_back(i);
					curv[i][2] = 1;
				}
			} else {
				if(curv[i][0] < curv[i-1][0]) {
					minima.push_back(i);
					curv[i][2] = 1;
				} else {
					curv[i][1] = 1;
				}
			}
		} else if(isClimbing) {
			if(curv[i][0] > curv[i+1][0]) {
				//is maxima
				//cout << "found maxima: " << curv[i-1][0] << endl;
				isClimbing = false;
				curv[i][1] = 1;
			}
		} else {
			if(curv[i][0] < curv[i+1][0]) {
				//is minima
				//cout << "found minima: " << curv[i-1][0] << endl;
				isClimbing = true;
				minima.push_back(i);
				curv[i][2] = 1;
			}
		}

	}

	//Sort minima array
	std::sort(std::begin(minima),
                std::end(minima),
                [&](int i1, int i2) { return curv[i1] < curv[i2]; } );
	cout << "MINIMA: " << curv[minima[0]][0] << endl;
	/*for(size_t k = 0; k < minima.size(); k++) {
		cout << curv[minima[k]][0] << " max: " << curv[minima[k]][1] << endl;
		cout << "min i: " << minima[k] << endl;

	}
	for(size_t j = 0; j < curv.size(); j++) {
		cout << curv[j][0] << " max: " << curv[j][1] << endl;

	}
	cout << endl;
	cout << endl;
	cout << "global max: " << curv[global_maxima][0] << " max: " << curv[global_maxima][1] << endl;
	cout << "global min: " << curv[global_minima][0] << " max: " << curv[global_minima][1] << endl;
	*/
}

/* Component besteht aus:
* - res[0] = Index des Minimums, dient gleichzeitig als Identifikator
* - res[1] = Index des Anfangspunktes
* - res[2] = Index des Endpunktes
* - res[3] = Index des Maximums
*/
template<class ElementType, class IndexType=vector<size_t>>
void calc_component(ElementType &curv, IndexType &res){

	//Keep adding smallest neighbours until maximum is found
	//cout << "Suche nach punkten in Minimum component: " << curv[res[0]][0] << endl;
	while(true)
	{

		if((res[1] >= 1) && (res[2] <= curv.size()-2)) {


			if(curv[res[1]-1][0] < curv[res[2]+1][0]) {
				res[1]--;
				//Maximum found left
				if(curv[res[1]][1] == 1) {
					res[3] = res[1];
					//cout <<  " Maximum left1: " << curv[res[1]][0] << endl;

					return;
				}
			} else {
				res[2]++;
				//Maximum found right
				if(curv[res[2]][1] == 1) {
					res[3] = res[2];
					//cout <<  " Maximum right1: " << curv[res[2]][0] << endl;
					return;
				}
			}
		} else if (res[1] >= 1) {
			res[1]--;
			//Maximum found left
			if((curv[res[1]][1] == 1) && (res[1] > 0)) {
				res[3] = res[1];
				cout <<  " Maximum left: " << curv[res[1]][0] << endl;

				return;
			}


		} else if (res[2] <= curv.size()-2) {
			res[2]++;
			//Maximum found right
			if((curv[res[2]][1] == 1) && (res[2]) < curv.size()-1) {

				res[3] = res[2];
				cout <<  " Maximum right: " << curv[res[2]][0] << " of minimum: " << curv[res[0]][0] << endl;

				return;
			}

		} else {
			return;
		}

	}
}

/*
* Berchnet die für den online algorithmus gebrauchten bars
* Problem zur Zeit noch: Jedes Mal wenn ein Maximum in zwei Komponenten gefunden wurde, setze ich das Maximum-Flag im
* curvature vector auf 0 füge das Minimum wieder in den minima vector ein. d.h. aber auch auch, dass die Komponente wieder von
* Anfang an aufgebaut wird.
* BarType = vector<vector<size_t>>
* ElementType = curvatures
*/
template <class ElementType, class IndexType=vector<size_t>, class BarType>
void calc_bars(ElementType &curv, IndexType minima, BarType &bars)
{
	vector<vector<size_t>> components;
	components.clear();
	bars.clear();
	//Gehe durch alle minima von geringstem bis höchstem
	for(size_t i = 0; i < minima.size(); i++){
		//Erstelle neue Komponente mit dem gefundenen Mimimumswert als Standardwert
		vector<size_t> new_comp(3);
		new_comp.clear();
		new_comp.push_back(minima[i]);
		new_comp.push_back(minima[i]);
		new_comp.push_back(minima[i]);
		new_comp.push_back(minima[i]);
		bool foundOld = false;
		bool foundMaxima = false;

		//Suche nach Minimum in bestehenden Komponenten. Setze derzeitige Komponente auf bestehende falls gefunden.
		for(size_t n = 0; n < components.size(); n++) {
			//cout << "Component: " << curv[components[n][0]][0] << " minimum: " << curv[minima[i]][0] << endl;
			if(components[n][0] == minima[i]) {
				//cout << "Found component with minimum" <<  curv[minima[i]][0] << endl;
				new_comp.assign(components[n].begin(), components[n].end());
				foundOld = true;
				components.erase(components.begin() + n);
			}
		}

		calc_component(curv, new_comp);
		//denotates if a component with same maximum was found
		//cout << "searching for components with same maximum..." << endl;
		for(size_t j = 0; j < components.size(); j++) {
			//cout << "New comp left: " << new_comp[1] << " right: " << new_comp[2] << endl;
			//cout << "Old comp left: " << components[j][1] << " right: " << components[j][2] << endl;
			if(new_comp[3] == components[j][3]) {
				//cout << "found same maximum: " << curv[new_comp[0]][0] << " and " << curv[components[j][0]][0] << endl;


				if(curv[new_comp[0]][0] < curv[components[j][0]][0]) {
					vector<size_t> new_bar = {components[j][0], components[j][3]};
					bars.push_back(new_bar);

					//Set left part of component
					if(components[j][0] > new_comp[0]) {
						components[j][1] = new_comp[1];
					//Set right part of component
					} else {
						components[j][2] = new_comp[2];
					}
					components[j][0] = new_comp[0];
					//minima.push_back(new_comp[0]);


				} else {
					vector<size_t> new_bar = {new_comp[0], new_comp[3]};
					bars.push_back(new_bar);

					//Set left part of component
					if(components[j][0] > new_comp[0]) {
						components[j][1] = new_comp[1];
					//Set right part of component
					} else {
						components[j][2] = new_comp[2];
					}
					//minima.push_back(components[j][0]);
				}
				minima.insert(minima.begin()+i+1, components[j][0]);

				//cout << "New Component is from " << curv[components[j][1]][0] << " to " << curv[components[j][2]][0] << endl;
				foundMaxima = true;
				//curv[components[j][0]][3] == 0;
				//the found components were merged, minimum has to be inserted into minima vector again

			}
		}
		if(!foundMaxima) {
			cout << "No component found, add new one: " << curv[new_comp[0]][0] << endl;
			components.push_back(new_comp);
		}
		if((i == minima.size()-1) && (components.size() > 1)) {
			cout << "last minimum reached, components left" << endl;
			//Set possible maximum Minimum to minimal Minimum
			size_t maxMin = minima[0];
			for(size_t n = 0; n < components.size(); n++) {
				cout << "component with minimum: " << curv[components[n][0]][0] << endl;
				if(curv[maxMin][0] < curv[components[n][0]][0]) {
					maxMin = components[n][0];
				}
			}
			cout << "insert min: " << curv[maxMin][0] << endl;
			minima.insert(minima.begin()+i+1, maxMin);
		}
		//cout << endl;

	}
	//Add last component as bar
	//bars.push_back({components[components.size()-1][0], components[components.size()-1][3]});
	for(size_t k = 0; k < components.size(); k++) {

		cout << "component #" << k << " start: " << curv[components[k][0]][0] << " end: " << curv[components[k][3]][0] << endl;
	}

	for (size_t f = 0; f < bars.size(); f++) {
		cout << "bar: "<< f << " left: " << curv[bars[f][0]][0] << " right: " << curv[bars[f][1]][0] << endl;
	}


}

//ONLY for testing purpuses
void initialiseTestVector(curvatures &test) {
	test.clear();
	test.push_back({0,0});
	test.push_back({-5,0});
	test.push_back({-8,0});
	test.push_back({15,0});
	test.push_back({30,0});
	test.push_back({-10,0});
	test.push_back({5,0});
	test.push_back({-40,0});
	test.push_back({-90,0});
	test.push_back({50,0});
	test.push_back({-12,0});
	test.push_back({13,0});
	test.push_back({-20,0});
	test.push_back({80,0});
	test.push_back({90,0});
	test.push_back({60,0});
	test.push_back({65,0});
	test.push_back({70,0});
	test.push_back({80,0});

}

template<class TrajectoryType, class ThresholdType=double>
TrajectoryType persistence_algo_off(TrajectoryType &traj, ThresholdType beta)
{
	//Schritt 1: Ableitung einer Kurvenfunktion


	//TODO: delete, only for testing
	curvatures curv(3);
	initialiseTestVector(curv);
	//Schritt 2: Berechne maxima und minima
	vector<size_t> minima;
	calc_extrema(curv, minima);
	cout << "done with extrema" << endl;
	for(size_t b = 0; b < curv.size(); b++) {
		cout << "Point: " << curv[b][0] << " min: " << curv[b][2] << " max: " << curv[b][1] << endl;
	}
	vector<vector<size_t> > bars;
	calc_bars(curv, minima, bars);

	bool startPoint = false;
	bool endPoint = false;
	vector<double> res(1);
	res.clear();
	//beta-persistent simplification
	cout << "bars size: " << bars.size() << endl;
	for(size_t i = 0; i < bars.size(); i++) {
		if((abs(curv[bars[i][0]][0]) + abs(curv[bars[i][1]][0])) < beta) {
			cout << "Beta smaller: " << abs(curv[bars[i][0]][0]) + abs(curv[bars[i][1]][0]) << endl;
			bars.erase(bars.begin() + i);
			i--;
		} else {
			cout << "Beta bigger" << endl;
			//Check if end points already are in simplification
			if(bars[i][0] == 0) {
				startPoint = true;
			}
			if(bars[i][1] == curv.size()-1) {
				cout << "end point found";
				endPoint == true;
			}
			res.push_back(bars[i][0]);
			res.push_back(bars[i][1]);

		}
	}
	sort(res.begin(), res.end());
	//cout << endl;
	//cout << "size: " << res.size() << endl;
	for(size_t j = 0; j < res.size(); j++) {
		cout << curv[res[j]][0] << endl;
	}

	if(!startPoint) {
		//Check points after start point and insert
		size_t toInsert = 1;
		bool insert = false;
		if(curv[res[0]][1] == 1){
			//first inserted extremum is maximum
			for(size_t k = 1; k < res[0]; k++) {
				if((curv[k][0] < curv[toInsert][0]) && (curv[k][2] == 1)) {
					//Inserted point between start and first inserted has to be minimum
					toInsert = k;
					insert = true;
				}
			}
		} else {
			//first insterted extremum is minimum
			for(size_t j = 1; j < res[0]; j++) {
				if((curv[j][0] > curv[toInsert][0]) && (curv[j][1] == 1)) {
					//Inserted point between start and first inserted has to be maximum
					toInsert = j;
					insert = true;
				}
			}
		}
		if(toInsert > -1) {
			res.insert(res.begin(), toInsert);
		}
		res.insert(res.begin(), 0);

	}
	if(!endPoint) {
		//Check points before an end point and insert
		size_t toInsert = curv.size()-2;
		bool insert = false;
		if(curv[res[res.size()-1]][1] == 1){
			//last inserted extremum is maximum
			for(size_t k = curv.size()-1; k > res[res.size()-1]; k--) {
				//cout << "check last points" << endl;
				if((curv[k][0] < curv[toInsert][0]) && (curv[k][1] == 1)) {
					toInsert = k;
					insert = true;
				}
			}
		} else {
			//last insterted extremum is minimum
			for(size_t j = curv.size()-1; j > res[res.size()-1]; j--) {
				if((curv[j][0] > curv[toInsert][0]) && (curv[j][2] == 1)) {
					toInsert = j;
					insert = true;
				}
			}
		}
		if(insert){
			res.insert(res.end(), toInsert);
		}
		res.insert(res.end(), curv.size()-1);
	}


	//res.insert(res.begin(), 0);
	//res.insert(res.end(), curv.size()-1);

	cout << endl << endl;

	for(size_t n = 0; n < res.size(); n++) {
		cout << curv[res[n]][0] << endl;
	}

	return traj;


}

//--------------------END ONLINE ALGORITHM-------------------------------

//--------------------START THRESHOLD SAMPLING---------------------------
/*
* Berechnet den Geschwindigkeitsvektor zwischen zwei Punkten
* returns vector<velocity, orientation>
*/
/*
vector<double> calc_velocity_vector(double lon1, double lon2, double lat1, double lat2, double t1, double t2)
{
	vector<double> velocity_vector;
	//velocity
	//cout << "P1: " << lon1 << "|" << lat1 << "|" << t1 << " - P2: " << lon2 << "|" << lat2 << "|" << t2 << endl;
	double mps = haversine(lat1, lat2, lon1, lon2) / ((t2 - t1) / 1000.0);
	//cout << "velocity: " << mps << endl;

	velocity_vector.push_back(mps);

	//orientation
	//double orientation = atan2((y2-y1),(x2-x1));

	//http://www.sunearthtools.com/tools/distance.php#txtDist_3
	double delta_lon = abs(lon1-lon2);
	if(delta_lon*(180/PI) > 180) {
		delta_lon = fmod(delta_lon*(180/PI), 180.0) *(PI/180);
	}
	double orientation = atan2(abs(lon1-lon2) ,log(tan((lat1/2) + (PI/4)) / tan((lat2/2) + (PI/4))));
	//cout << "orientation: " << orientation*(180/PI) << endl;

	velocity_vector.push_back(orientation);
	return velocity_vector;
}


template<class TrajectoryType, class TimesType>
TrajectoryType threshold_sampling(TrajectoryType &traj, TimesType &times, double velocity_thresh, double orientation_thresh)
{
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
		cout << "result size: " << resultSize << endl;

		//START Calculate upper and lower angle and radius for trajectory and sample based last two points.

		//For velocity vector of trajectory, take last 2 elements in original trajectory
		vector<double> velocity_vector_traj = calc_velocity_vector(traj[i-2][0], traj[i-1][0],
											traj[i-2][1], traj[i-1][1], times[i-2], times[i-1]);
		//For velocity vector of sample, take last 2 elements in result trajectory
		vector<double> velocity_vector_sample = calc_velocity_vector(result[resultSize-2][0], result[resultSize-1][0], result[resultSize-2][1],
										result[resultSize-1][1], resultTimes[resultSize-2], resultTimes[resultSize-1]);


		//Calculate the radius for trajectory and sample based using the calculated velocity vectors
		double upper_radius_traj = ((times[i] - times[i-1]) / 1000.0) *
										(velocity_vector_traj[0]*(1 + velocity_thresh));
		double lower_radius_traj = ((times[i] - times[i-1]) / 1000.0) *
										(velocity_vector_traj[0]*(1 - velocity_thresh));
		cout << "upper radius_traj : " << upper_radius_traj << endl;
		cout << "lower radius_traj : " << lower_radius_traj << endl;

		double upper_radius_sample = ((times[i] - times[i-1]) / 1000.0) *
										(velocity_vector_sample[0]*(1 + velocity_thresh));
		double lower_radius_sample = ((times[i] - times[i-1]) / 1000.0) *
										(velocity_vector_sample[0]*(1 - velocity_thresh));
		cout << "upper_radius_sample : " << upper_radius_sample << endl;
		cout << "lower_radius_sample : " << lower_radius_sample << endl;

		//Calculate circle planes for trajectory and sample
		//double circle_plane_traj = (pow(upper_radius_traj, 2) * PI) - (pow(lower_radius_traj, 2)*PI);
		//cout << "circle_plane_traj: " << circle_plane_traj << endl;

		//double circle_plane_sample = (pow(upper_radius_sample, 2) * PI) - (pow(lower_radius_sample, 2)*PI);
		//cout << "circle_plane_sample: " << circle_plane_sample << endl;

		//Calculate orientation angle for both from given slopes
		double upper_slope_traj = velocity_vector_traj[1]*(180/PI)  + orientation_thresh;
		double lower_slope_traj = velocity_vector_traj[1]*(180/PI)  - orientation_thresh;
		//double angle_traj = atan(abs((upper_slope_traj - lower_slope_traj) /
		//						(1 + upper_slope_traj*lower_slope_traj))) * (180/PI);
		cout << "upper angle traj: " << upper_slope_traj << endl;
		cout << "lower angle traj: " << lower_slope_traj << endl;

		double upper_slope_sample = velocity_vector_sample[1]*(180/PI)  + orientation_thresh;
		double lower_slope_sample = velocity_vector_sample[1]*(180/PI)  - orientation_thresh;
		//double angle_sample = atan(abs((upper_slope_sample - lower_slope_sample) /
		//						(1 + upper_slope_sample*lower_slope_sample))) * (180/PI);
		cout << "upper_slope_sample: " << upper_slope_sample << endl;
		cout << "lower_slope_sample: " << lower_slope_sample << endl;

		//END

		//START Calculate angle and radius for point to insert
		vector<double> velocity_vector_insert = calc_velocity_vector(traj[i-1][0], traj[i][0],
											traj[i-1][1], traj[i][1], times[i-1], times[i]);

		double radius_insert = ((times[i] - times[i-1]) / 1000.0) *
										velocity_vector_insert[0];
		cout << "radius_insert: " << radius_insert << endl;

		double angle_insert = velocity_vector_insert[1] * (180/PI);
		//END

		//Check if new points angle and radius is inside of the calculated
		cout << "angle insert: " << angle_insert << endl;
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



		//double safe_area_traj = (angle_traj / 360) * circle_plane_traj;
		//cout << "safe area: " << safe_area_traj << endl;
		//double safe_area_sample = (angle_sample / 360) * circle_plane_sample;

		//TODO: Wie kann ich bestimmen, ob ein Punkt in der Ebene liegt?
		//Gegeben:
		//	- Ausgangspunkt: traj[i]
		//  - Fläche der Ebene: safe_area_traj
		//	- Richtung (in Form von Winkel): angle_traj
		//	- Geschwindigkeit: velocity_vector[0]
		// Idee:
		//	- Berechne Geschwindigkeitsvektor von vorherigen und neuem Punkt
		// 	- Berechne Winkel und die beiden radien ()
	}

	return result;


}
*/

//TODO Implement STTrace
template<class trajectory>
trajectory STTrace(trajectory &traj, trajectory &times, double velocity_thresh, double orientation_thresh)
{

}

//----------------------END THRESHOLD SAMPLING----------------------------------------

//----------------------START DBSCAN DENSITY ALGO-------------------------------------


/*
* Used for executing a command in bash and returning the output
*/
std::string exec(const char* cmd) {
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


/*
* Density based segmentation using DBSCAN
* Gets loaded trajectory and times in order to sample them first (will be saved in folder /sample of chosen user folder)
*
* Proceeds to call ELKI DBSCAN implementation with epsilon and minPts parameters given by user, saves resulting clusters 7
* in /cluster folder of chosen user
*
* Segments trajectory according to resulting clusters where the most middle point is chosen for segmentation
*/
template<class TrajectoryType>
size_t density_dbscan_segmentation(TrajectoryType &traj, TrajectoryType &res, double epsilon, int minPts, const std::string &dataName, const std::string &pathToData)
{
	//Schritt 1: sample data (Evtl auch davor in der main aufrufen)

	//Schritt 2: Rufe DBSCAN auf übergebenen Datenpfad auf
	std::stringstream cmd;
	cmd << "java -jar ../../elki.jar KDDCLIApplication -dbc.in " << pathToData << " -algorithm clustering.DBSCAN " <<
	"-dbscan.epsilon " << epsilon << " -dbscan.minpts " << minPts << " -resulthandler ResultWriter -out clustering/" << dataName;
	string s = cmd.str();
	const char* p = s.c_str();
	string result = exec(p);

	cout << result << endl;

	//Erster Punkt von Trajektorie muss im Ergebnis sein
	res.clear();
	res.push_back({traj[0][0], traj[0][1]});

	vector<vector<double> > coords(3);


	//Gehe gefundene Cluster in Ordner durch
	int i = 0;
	while(true)
	{
		stringstream path_name;
		path_name << "clustering/" << dataName << "/Cluster_" << i << ".txt";
		ifstream file(path_name.str().c_str());
		if(!file) {
			cout << "No more clusters found... " << endl;
			break;
		}
		cout << "Found file " << endl;
		coords.clear();
		string line;
		while(getline(file, line)) {
			size_t found = line.find("#");

			if(found != string::npos){
				cout << "Found commentary... " << endl;
				continue;
			} else {
				cout << "Found clusterline" << endl;
			}
			istringstream iss(line);
			string IDString;
			double lat, lon;
			if(iss >> IDString >> lat >> lon) {
			}
			IDString = IDString.substr(IDString.find("=")+1, string::npos);
			cout << "ID " << IDString << endl;
			iss.str(IDString);
			double ID;
			iss >> ID;

			coords.push_back({ID,lat,lon});
		}
		cout << "cluster: " << coords[0][0] << " " << coords[0][1] << " " <<  coords[0][2] << endl;
		//Schritt 3: Suche repräsentierbare Punkte in Cluster
		//Theoretisch könnte man hier noch den Punkt als Segmentierung wählen, der am nächsten in der MItte ist.
		res.push_back(traj[coords[0][0]]);

    	file.close();
		i++;
	}

	//Letzter Punkt in Trajektorie muss im Ergebnis sein
	res.push_back(traj[traj.size()-1]);


	return res.size();

}


//----------------------END DBSCAN DENSITY ALGO---------------------------------------

double latitude_converter(std::string lat)
{
	int length = lat.length();
	std::stringstream formated_lat;
	switch (length)
	{
		case 8:
			formated_lat << lat.substr(0,1) << "." << lat.substr(1,8);
			break;
		case 9:
			formated_lat << lat.substr(0,2) << "." << lat.substr(2,8);
			break;
		default:
			return 0.0;
	}
	double result;
	formated_lat >> result;
	return result;
}

double longitude_converter(std::string lon)
{
	int length = lon.length();
	std::stringstream formated_lon;
	switch (length)
	{
		case 8:
			formated_lon << lon.substr(0,1) << "." << lon.substr(1,8);
			break;
		case 9:
			formated_lon << lon.substr(0,2) << "." << lon.substr(2,8);
			break;
		case 10:
			formated_lon << lon.substr(0,3) << "." << lon.substr(3,8);
			break;
		default:
			return 0.0;
	}
	double result;
	formated_lon >> result;
	return result;
}

void export_google_times(std::string &data_path, std::string &user_name)
{
	std::cout << std::fixed;
	try
       {

		   time_t rawtime;
		   time (&rawtime);
		   struct tm * lastTime;
		   lastTime = localtime(&rawtime);
			char lastBuffer[256];
			strftime(lastBuffer, sizeof(lastBuffer), "%d", lastTime);
			cout << "Current day: " << std::string(lastBuffer) << endl;

			std::stringstream data_location;
			data_location << "data/processed/GoogleNow/" << user_name << "/";
			std::stringstream times_location;
			times_location << "data/processed/GoogleNow/" << user_name << "/times/";
			string data_save_path = data_location.str();
			string time_save_path = times_location.str();
			size_t i = 0;
			data_location.str("");
			times_location.str("");

			boost::property_tree::ptree pt;
			read_json(data_path, pt);
           BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child("locations"))
           {
			   std::string dateName = std::string(lastBuffer);
			   std::string curTime = "";
			   double lat = 0.0;
			   double lon = 0.0;
			   BOOST_FOREACH(boost::property_tree::ptree::value_type &p, v.second)
               {
				   if(p.first == "timestampMs") {
					   curTime = p.second.get_value<std::string>().substr(0,10);
					   std::stringstream tmp(curTime);
					   long rawtime;
					   tmp >> rawtime;
					   const time_t tmptime = (const time_t) rawtime;
					   struct tm * curTime;
					   curTime = localtime(&tmptime);

					   char curBuffer[256];
					   strftime(curBuffer, sizeof(curBuffer), "%d", curTime);
					   if(std::string(curBuffer) != std::string(lastBuffer)) {
						   strcpy(lastBuffer, curBuffer);
						   dateName = std::string(curBuffer);
						   //lastBuffer = curBuffer;
						   cout << "Found new date" << endl;
						   i++;
					   }
				   } else if (p.first == "latitudeE7") {
					   lat = latitude_converter(p.second.get_value<std::string>());
				   } else if (p.first == "longitudeE7") {
					   lon = longitude_converter(p.second.get_value<std::string>());
				   }
			   }
			   //Write to files
			data_location << data_save_path  << i << ".csv";
			times_location << time_save_path << i << ".csv";

			std::ofstream traj_out;
			traj_out << fixed;
			traj_out.precision(7);
			traj_out.open(data_location.str().c_str(), ios::app);

   			std::ofstream time_out;
			time_out.open(times_location.str().c_str(), ios::app);
   			time_out.precision(13);
   			if (!traj_out.is_open() || !time_out.is_open())
   			{
   				throw(std::runtime_error("Could not open output file."));
   			}
			traj_out << lat << " " << lon << std::endl;
			time_out << curTime << std::endl;
			data_location.str("");
			times_location.str("");
			traj_out.close();
			time_out.close();
           }

       }
       catch (std::exception const& e)
       {
           std::cerr << e.what() << std::endl;
       }
}

int main(int argc, char *argv[])
{
	std::cout << std::setprecision(10);
	cout << "starting" << endl;
	string dataPath = "data/raw/testing.json";
	string username = "testing";
	export_google_times(dataPath, username);
	   return 0;

	trajectory traj(2);
	vector<double> times(1);
	trajectory resample_result(2);
	vector<double> timesRes(1);

	string pathToData = "../../python-server-sample/TestTraj.csv";
	loadXY(pathToData, traj);
	loadTime("../../python-server-sample/TestTrajTimes.csv", times);
	std::cout << "loaded data" << endl;

	trajcomp::tools::resample_trajectory(traj, times, resample_result,
		timesRes, "/home/jan/Desktop/Masterarbeit/trajcomputing_code/v_2/libtrajcomp-src/doc/md/resample/", "TestTraj", 2.0);
	cout << "original length: " << traj.size() <<  " resample: " << resample_result.size() << endl;
	std::cout << "resampled" << endl;

	return 0;

	trajectory input(2);
	string dataName = "Test";
	//string pathToData = "../../python-server-sample/TestTraj.csv";
	string path_to_elki = "/home/jan/elki.jar";
	loadXY("../../python-server-sample/TestTraj.csv", input);

	trajectory dbscan_based = trajcomp::DBSCAN_segmentation(input, 0.001,
		2, dataName, pathToData, path_to_elki);

	//density_dbscan_segmentation(input, dbscan_based, 0.00000002, 1, dataName, pathToData);
	cout << "Result of dbscanbased: " << dbscan_based.size() << endl;
	for (size_t j = 0; j < dbscan_based.size(); j++) {
		cout << dbscan_based[j][0] << " " << dbscan_based[j][1] << endl;
	}
	return 0;


	trajectory online(2);

	persistence_algo_off(online, 25.0);
	return 0;

	trajectory test(2);
	vector<double> ti(1);
	trajectory threshold(2);
	//Muss im Format - {longitute latitude} - sein
	loadXY("../../python-server-sample/TestTraj.csv", test);
	//Muss im Format - {Millisekunden} - sein
	loadTime("../../python-server-sample/TestTrajTimes.csv", ti);
	threshold = trajcomp::threshold_sampling(test, ti, 0.7, 20.0);
	cout << "original size: " << test.size() << endl;

	cout << "threshold size: " << threshold.size() << endl;

	return 0;

	trajectory dp_traj(2);
	trajectory googleNowTraj(2);
	loadXY("../../python-server-sample/TestTraj.csv", googleNowTraj);
	//loadXY("./data/processed/GoogleNow/User1/trajectory3.dat", googleNowTraj);

	dp_traj = trajcomp::douglas_peucker(googleNowTraj, 0.001);
//	dp_traj.summary();
	cout << "original size: " << googleNowTraj.size() << endl;

	cout << "DP size: " << dp_traj.size() << endl;
	//timesTraj.load("./data/processed/GoogleNow/User1/times3.dat");
	//timesTraj.summary();
	//cout << "Density: " << endl;
	//trajectory segmentedTraj = density_segmentation(googleNowTraj, 1.0);
	//segmentedTraj.dump();
	//cout << "Douglas Peucker" << endl;
	//cout << trajcomp::tools::make_string2(segmentedTraj);

	//cout << "d(test,test2) = "<< trajcomp::discrete_frechet(test,test2)<< endl;

	//cout << "d(t2,t2) = "<< trajcomp::discrete_frechet(googleNowTraj,dp_traj)<< endl;
	return 0;
}
