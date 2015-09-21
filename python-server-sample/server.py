#from methods import *
import bottle
from bottle import request, route, hook, response
import json
import backendclient;
import beaker.middleware
import methods

session_opts = {
    'session.type': 'file',
    'session.data_dir': './session/',
    'session.auto': True,
}

app = beaker.middleware.SessionMiddleware(bottle.app(), session_opts)

@hook('before_request')
def setup_request():
	request.session = request.environ['beaker.session']
	try:
		if 'messages' in request.session:
			request.messages = json.loads(request.session['messages']);
		else:
			request.messages=[];
	except ValueError:
		request.messages=[];



#do_restore();

@route('/')
def base():
	response.headers['Content-Type'] = 'text/html'
	f = open('template.html', 'r')
	return f.read()

@route('/distance')
def get_distance():
	print "done some clicking"
	return "{}"

@route('/click')
def server_clicked():
	x = float(request.query.lat)
	y = float(request.query.lon)
	print "Click"+str(x) + " " + str(y)
	request.messages.append(['click',x,y]);
	print request.messages;
	request.session['messages']=json.dumps(request.messages);
	return "{}"

@route('/polylines')
def get_polylines():
	# Format is array of arrays
	p1 = [[11.55800670,48.14483530], [12.55800670,48.14483530],[12.55800670,49.14483530],[11.55800670,49.14483530]];
	p2 = [[13.55800670,48.14483530], [14.55800670,48.14483530],[13.55800670,49.14483530],[13.55800670,49.14483530]];

	return json.dumps([p1,p2]);

@route('/view')
def view():
	response.headers['Content-Type'] = 'image/png'
	return methods.handle_events(request);

@route('/get/<theid>')
def get_by_id(theid='0'):
	return methods.get_trajectory(theid);

@route('/intersect')
def server_clicked():
	intersects = json.loads(request.query.hashes)
	return methods.handle_intersects(intersects);

@route('/jaccard')
def jaccard_clicked():
	intersects = json.loads(request.query.hashes)
	return methods.handle_jaccard(intersects);

@route('/subsetcorr')
def subsetcorr_clicked():
	intersects = json.loads(request.query.hashes)
	a = backendclient.backend_subsetcorr(intersects);
	print a
	return json.dumps(a);

@route('/numjaccard/<theid>')
def jaccard_id(theid='0'):
	return methods.get_jaccard_by_id(theid);


@route('/reset')
def reset():
	request.session.delete()
	return "{}"

@route('/douglasPeucker')
def save():
    epsilon = request.query.eps;
    coords = json.loads(request.query.coords);
    target = open('TestTraj.csv', 'w');
    target.truncate();
    for key in coords:
        target.write(str(key));
        target.write("\n");

    target.close();
    # Send location of save trajectory to backend server
    return methods.get_douglas_peucker("../TestTraj.csv", epsilon);
	#return methods.do_DBSCAN("TestTraj.csv");

@route('/dbscan')
def dbscan():
    epsilon = request.query.eps;
    minPts = request.query.minPts;
    coords = json.loads(request.query.coords);
    target = open('TestTraj.csv', 'w');
    target.truncate();
    for key in coords:
        target.write(str(key));
        target.write("\n");

    target.close();
    # Send location of save trajectory to backend server
    return methods.get_DBSCAN("../TestTraj.csv", minPts, epsilon);


@route('/threshold')
def dbscan():
    vel_thresh = request.query.vel;
    or_thresh = request.query.orientation;
    coords = json.loads(request.query.coords);
    times = json.loads(request.query.times);

    target = open('TestTraj.csv', 'w');
    target.truncate();
    for key in coords:
        target.write(str(key));
        target.write("\n");
    target.close();

    target = open('TestTrajTimes.csv', 'w');
    target.truncate();
    for time in times:
        target.write(str(time));
        target.write("\n");
    target.close();
    return methods.get_threshold_data(vel_thresh, or_thresh);

@route('/resampleData')
def resample():
    coords = json.loads(request.query.coords);
    dist_for_resampling = request.query.dist;
    times = json.loads(request.query.times);
    target = open('TestTraj.csv', 'w');
    target.truncate();
    for key in coords:
        target.write(key);
        target.write("\n");
    target.close();
    target = open('TestTrajTimes.csv', 'w');
    target.truncate();
    for time in times:
        target.write(str(time));
        target.write("\n");
    target.close();
    # Send location of save trajectory to backend server
    R = methods.resample_data("TestTraj", dist_for_resampling);
    print R;
    return R;

#@get('/query/<thefilter>')
#def query_filter( thefilter="NoFilter"):
#    return do_query_filter(thefilter);

#@put('/store')
#def store_filter( ):
#    return do_put_filter(json.load(request.body));


#@get('/clear')
#def query_clear( ):
#    return do_clear();


#run(host='trajcomp.mobile.ifi.lmu.de', port=8080)
bottle.run (app=app,host='localhost',port=8081);

print "Started server"
#print "Saving state:"
#do_save();
