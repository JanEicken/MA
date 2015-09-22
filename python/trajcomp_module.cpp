#include<Python.h>
#include<stdlib.h>
#include<trajcomp/trajcomp_all.hpp>

#include<vector>
#include<string>

// PyExc_IOError
// PyExc_RuntimeError

/*Essentially, we have a store of trajectories
 * held in the library here
 * */

#define GLOBAL_BOUND_CHECK(x) { if (x < 0 || x >= trajectories.size()) \
			{ \
				PyErr_SetString(PyExc_IndexError,"Trajcomp: Out of bounds");\
				return NULL; \
			} }



typedef trajcomp::trajectory<double> trajectory;
typedef vector<double> times;
trajcomp::default_segment_distance<trajectory::value_type> dseg;


std::vector<trajectory> trajectories;
std::vector<trajectory> google_now_trajectories;
std::vector<trajectory> google_now_times;

const std::string resample_location = "../resample/";
const std::string resample_name = "TestTraj";


// load a trajectory
static PyObject *
trajcomp_load(PyObject *self, PyObject *args)
{
	const char *argument;
	if (!PyArg_ParseTuple(args, "s", &argument))
        return NULL;
    trajectory t;
    std::cout << "loading" << argument << std::endl;
    try{
    t.load(std::string(argument));
   }catch(std::runtime_error e)
   {
	   PyErr_SetString(PyExc_IOError,"Could not read file");
	   return NULL;
   }

    trajectories.push_back(t);
    return Py_BuildValue("i", trajectories.size() -1);
}

static PyObject *
trajcomp_summary(PyObject *self, PyObject *args)
{
	unsigned int i;
	if (!PyArg_ParseTuple(args, "i", &i))
        return NULL;
    GLOBAL_BOUND_CHECK(i);
    trajectories[i].summary();

	return Py_BuildValue("i", trajectories[i].size());
}

static PyObject *
trajcomp_get(PyObject *self, PyObject *args)
{
	unsigned int i;
	if (!PyArg_ParseTuple(args, "i", &i))
        return NULL;
    GLOBAL_BOUND_CHECK(i);

    PyObject *list = PyList_New(0);
    trajectory::Base::iterator it;
    for (it = trajectories[i].begin();
		 it != trajectories[i].end();
		  it++)
    {
		PyObject *child = PyList_New(0);
		trajectory::ElementType::iterator it2;
		for (it2 = (*it).begin();
			it2 != (*it).end();
			it2++)
		{
				PyList_Append(child, Py_BuildValue("d",*it2));
		}
		PyList_Append(list,child);
    }
	return list;
}


static PyObject *
trajcomp_uniform_select(PyObject *self, PyObject *args)
{
	unsigned int i;
	unsigned int k;
	unsigned int last;

	if (!PyArg_ParseTuple(args, "iii", &i, &k, &last))
        return NULL;
    GLOBAL_BOUND_CHECK(i);

    trajectory t = trajcomp::uniform_select(trajectories[i],k,(last == 1));
    trajectories.push_back(t);
    return Py_BuildValue("i", trajectories.size() -1);
}


static PyObject *
trajcomp_douglas_peucker(PyObject *self, PyObject *args)
{
	unsigned int i;
	double epsilon;
	if (!PyArg_ParseTuple(args, "id", &i, &epsilon))
        return NULL;
    GLOBAL_BOUND_CHECK(i);
    trajectory t = trajcomp::douglas_peucker(trajectories[i],epsilon,dseg);
    trajectories.push_back(t);
    return Py_BuildValue("i", trajectories.size() -1);
}



/*Frechet Distance*/

static PyObject *
trajcomp_frechet(PyObject *self, PyObject *args)
{
	unsigned int i,j;
	double epsilon;
	if (!PyArg_ParseTuple(args, "iid", &i, &j, &epsilon))
        return NULL;
    GLOBAL_BOUND_CHECK(i);
    GLOBAL_BOUND_CHECK(j);
    double ret = trajcomp::frechet(trajectories[i],trajectories[j],epsilon);
    return Py_BuildValue("d", ret);
}


static PyObject *
trajcomp_point_segment(PyObject *self, PyObject *args)
{
	double x,y,x2,y2,x3,y3;

	if (!PyArg_ParseTuple(args, "dddddd", &x,&y,&x2,&y2,&x3,&y3))
        return NULL;
    std::vector<double> u; u.push_back(x);  u.push_back(y);
    std::vector<double> v; v.push_back(x2); v.push_back(y2);
    std::vector<double> p; p.push_back(x3); p.push_back(y3);

    trajcomp::default_segment_distance<std::vector<double>> d;

    double ret = d(u,v,p);
    return Py_BuildValue("d", ret);
}



static PyObject *
trajcomp_geolife(PyObject *self, PyObject *args)
{
	unsigned int i;
	const char *s;
	if (!PyArg_ParseTuple(args, "is", &i,&s))
        return NULL;
    // Load the first i trajectories into the list and return the first handle

    int first_index = trajectories.size();

    trajcomp::geolife_reader<trajectory> reader(s,i);
	std::vector<trajectory> data;
	std::cout << "Directories enumerated: " << reader.load(data) << std::endl;

    trajectories.insert(trajectories.end(),data.begin(),data.end());
    std::cout << "Trajectories now " << trajectories.size() << std::endl;
	std::cout << "First trajectory: " << trajectories[0][0][0] << std::endl;
	for(size_t i=0; i < trajectories[0][0].size(); i++) {
		std::cout << trajectories[0][0][i] << " | ";
	}
	std::cout << std::endl << "blubb" << std::endl;
    return Py_BuildValue("i", first_index);
}


/*
static PyObject *
trajcomp_googleNow(PyObject *self, PyObject *args)
{
	unsigned int i;
	const char *s;
	if (!PyArg_ParseTuple(args, "is", &i, &s)){
		std::cout << "Failure in py args" << endl;
		return NULL;
	}

	int first_index = trajectories.size();

    trajcomp::googleNow_reader<trajectory> reader(s,i);
	std::vector<trajectory> data;
	std::cout << "Directories enumerated: " << reader.load(data, false) << std::endl;

    google_now_trajectories.insert(trajectories.end(),data.begin(),data.end());
    std::cout << "Trajectories now " << google_now_trajectories.size() << std::endl;
	std::cout << "First trajectory: " << google_now_trajectories[0][0][0] << std::endl;

	std::vector<times> dataTimes;
	reader.load(dataTimes, true);
	std::cout << std::endl << "blubb" << std::endl;
    return Py_BuildValue("i", first_index);


}
*/

static PyObject *
trajcomp_size(PyObject *self, PyObject *args)
{
	return Py_BuildValue("i", trajectories.size());
}


//TODO: Function to use trajcomp threshold guided sampling
// Input: Path to trajectory
static PyObject *
trajcomp_threshold_sampling(PyObject *self, PyObject *args)
{
  const char *s;
	if(!PyArg_ParseTuple(args, "s", &s))
			return NULL;

	return Py_BuildValue("i", trajectories.size());
}

// Input: Path to trajectory, epsilon
static PyObject *
trajcomp_douglas_peucker_online(PyObject *self, PyObject *args)
{
	const char *s;
	//const char *epsilon;
	double epsilon;
	if (!PyArg_ParseTuple(args, "sd", &s, &epsilon)) {
		std::cout << "Error in parsing params from python to c++" << std::endl;
		return NULL;
	}
	trajectory traj(2);
	trajectory dp_result(2);
	string data(s);
//	string data2(epsilon);
 	loadXY(data, traj);
	std::cout << "Douglas Peucker original length: " << traj.size() << endl;
	std::cout << "Epsilon: " << epsilon << endl;
 	dp_result = trajcomp::douglas_peucker(traj, epsilon);
	std::cout << "Douglas peucker  result length: " << dp_result.size() << endl;;
 	PyObject *list = PyList_New(0);
 	trajectory::Base::iterator it;
 	for (it = dp_result.begin(); it != dp_result.end(); it++)
	{

		PyObject *child = PyList_New(0);
		trajectory::ElementType::iterator it2;
		for (it2 = (*it).begin(); it2 != (*it).end(); it2++)
		{
			std::cout << *it2;
			PyList_Append(child, Py_BuildValue("d",*it2));
		}
		std::cout << std::endl;
		PyList_Append(list,child);
	}
	return list;
}

static PyObject *
trajcomp_DBSCAN_segmentation(PyObject *self, PyObject *args)
{
	std::cout << "enter python module dbscan" << endl;
	//Name of the folder in Which the clusters will be saved
	const char *d;
	const char *n;
	const char *p;
	double epsilon;
	unsigned int minPts;

	if (!PyArg_ParseTuple(args, "disss", &epsilon, &minPts, &n, &d, &p)) {
		std::cout << "Error in parsing params from python to c++" << std::endl;
		return NULL;
	}
	string dataPath(d);
	string pathToElki(p);
	// e.g. /googleNow/Jan/1/
	string clusterName(n);
	trajectory traj(2);
	loadXY(dataPath, traj);
	trajectory dbscan_result = trajcomp::DBSCAN_segmentation(traj, epsilon, minPts, clusterName, dataPath, pathToElki);

	PyObject *list = PyList_New(0);
 	trajectory::Base::iterator it;
 	for (it = dbscan_result.begin(); it != dbscan_result.end(); it++)
	{

		PyObject *child = PyList_New(0);
		trajectory::ElementType::iterator it2;
		for (it2 = (*it).begin(); it2 != (*it).end(); it2++)
		{
			//std::cout << *it2 << " ";
			PyList_Append(child, Py_BuildValue("d",*it2));
		}
		//std::cout << std::endl;
		PyList_Append(list,child);
	}
	return list;
}

static PyObject *
trajcomp_resample(PyObject *self, PyObject *args)
{
	std::cout << "enter python module resample" << endl;
	const char *d;
	double h;
	const char *p;
	if (!PyArg_ParseTuple(args, "sds", &d, &h, &p)) {
		std::cout << "Error in parsing params from python to c++" << std::endl;
		return NULL;
	}
	string dataName(d);
	string dataPath(p);
	trajectory traj(2);
	vector<double> times(1);
	trajectory resample_result(2);
	vector<double> timesRes(1);
	loadXY("../TestTraj.csv", traj);
	loadTime("../TestTrajTimes.csv", times);
	std::cout << "loaded data" << endl;
	std::cout << "Data name and path: " << dataName << " | " << dataPath << endl;
	trajcomp::tools::resample_trajectory(traj, times, resample_result,
		timesRes, dataPath, dataName, h);

	std::cout << "resampled trajectory" << endl;
	PyObject *traj_list = PyList_New(0);
	PyObject *list = PyList_New(0);
 	trajectory::Base::iterator it;
 	for (it = resample_result.begin(); it != resample_result.end(); it++)
	{

		PyObject *child = PyList_New(0);
		trajectory::ElementType::iterator it2;
		for (it2 = (*it).begin(); it2 != (*it).end(); it2++)
		{
			std::cout << *it2 << " ";
			PyList_Append(child, Py_BuildValue("d",*it2));
		}
		PyList_Append(traj_list,child);
	}
	PyList_Append(list, traj_list);

	std::vector<double>::iterator it3;
	PyObject *times_list = PyList_New(0);
	for (it3 = timesRes.begin(); it3 != timesRes.end(); it3++)
	{

		PyList_Append(times_list, Py_BuildValue("d",*it3));
	}
	PyList_Append(list, times_list);
	return list;
}

static PyObject *
trajcomp_threshold(PyObject *self, PyObject *args)
{
	std::cout << "enter python module threshold" << endl;
	const char* d;
	const char* t;
	double v;
	double o;
	if (!PyArg_ParseTuple(args, "ddss", &v, &o, &d, &t)) {
		std::cout << "Error in parsing params from python to c++" << std::endl;
		return NULL;
	}
	trajectory traj(2);
	trajectory result(2);
	std::string dataPath(d);
	std::string timesPath(t);
	vector<double> times(1);
	trajectory resample_result(2);
	loadXY(dataPath, traj);
	loadTime(timesPath, times);
	result = trajcomp::threshold_sampling(traj, times, v, o);

	PyObject *list = PyList_New(0);
	trajectory::Base::iterator it;
 	for (it = result.begin(); it != result.end(); it++)
	{

		PyObject *child = PyList_New(0);
		trajectory::ElementType::iterator it2;
		for (it2 = (*it).begin(); it2 != (*it).end(); it2++)
		{
			std::cout << *it2 << " ";
			PyList_Append(child, Py_BuildValue("d",*it2));
		}
		std::cout << std::endl;
		PyList_Append(list,child);
	}
	return list;
}


static PyObject *
trajcomp_system(PyObject *self, PyObject *args);


static PyMethodDef TrajcompMethods[] = {
    {"system",  trajcomp_system, METH_VARARGS, "Execute a shell command."},
    {"load",  		trajcomp_load, METH_VARARGS, "Load a trajectory."},
    {"summary",  	trajcomp_summary, METH_VARARGS, "Summary of a trajectory."},
    {"get",  	trajcomp_get, METH_VARARGS, "Return trajectory as list."},
    {"uniform_select",  	trajcomp_uniform_select, METH_VARARGS, "Simplify using uniform select."},
    {"douglas_peucker",  	trajcomp_douglas_peucker, METH_VARARGS, "Simplify using Douglas Peucker."},
	{"douglas_peucker_online",  trajcomp_douglas_peucker_online, METH_VARARGS, "Simplify using Douglas Peucker online."},
    {"frechet",  	trajcomp_frechet, METH_VARARGS, "Calculate Frechet Distance (eps <0) or solve decision problem (eps>0)."},
    {"geolife",  	trajcomp_geolife, METH_VARARGS, "Load Geolife Data."},
    {"pld",  	trajcomp_point_segment, METH_VARARGS, "Point Segment Distance [debug]."},
    {"size",  	trajcomp_size, METH_VARARGS, "Count trajectories."},
	{"threshold_sampling", trajcomp_threshold_sampling, METH_VARARGS, "Use threshold guided sampling on input trajectory"},
	{"dbscan", trajcomp_DBSCAN_segmentation, METH_VARARGS, "DBSCAN segmentation on input trajectory path"},
	{"resample", trajcomp_resample, METH_VARARGS, "Resample input trajectory"},
	{"threshold", trajcomp_threshold, METH_VARARGS, "Calculate threshold sampling trajectory"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyObject *
trajcomp_system(PyObject *self, PyObject *args)
{
    const char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return Py_BuildValue("i", sts);
}

PyMODINIT_FUNC
inittrajcomp(void)
{
    PyObject *m;

    m = Py_InitModule("trajcomp", TrajcompMethods);
    if (m == NULL)
		std::cout << "fault" << endl;
        return ;

    //SpamError = PyErr_NewException("spam.error", NULL, NULL);
    //Py_INCREF(SpamError);
    //PyModule_AddObject(m, "error", SpamError);
}
