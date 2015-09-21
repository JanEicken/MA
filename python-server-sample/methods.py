#import trajcomp;
from operator import itemgetter
#import matplotlib.pyplot as plt
#from matplotlib.patches import Rectangle
#from base64 import *;
#import tempfile;
#import Image, ImageDraw;

import backendclient;
import json;
import mercator;

def get_trajectory(theid):
	a = backendclient.backend_get(theid);
	return json.dumps(a);


def handle_intersects(hashes):
	a = backendclient.backend_intersect(hashes);
	return json.dumps(a);

def handle_jaccard(hashes):
	a = backendclient.backend_jaccard(hashes);
	return json.dumps(a);

def get_jaccard_by_id(theid):
	a = backendclient.backend_getjaccard(theid);
	return json.dumps(a);


def get_DBSCAN(dataName, minPts, eps):
	a = backendclient.backend_DBSCAN_segmentation(dataName, minPts, eps);
	return json.dumps(a);

def get_douglas_peucker(dataName, eps):
	print "backend eps";
	print eps;
	a = backendclient.backend_DP_segmentation(dataName, eps);
	print a;
	return json.dumps(a);

def resample_data(dataName, dist):
	print dataName;
	print dist;
	a = backendclient.backend_resample_data(dataName, dist);
	return json.dumps(a);

def get_threshold_data(vel, orientation):
	a = backendclient.backend_threshold_data(vel, orientation);
	return json.dumps(a);


cfg_width=500
cfg_height=500





# RemoteCanvas Implementation

def handle_events(request):
	evt = request.messages;
	for e in evt:
		if e[0] == "click":
			handle_click(request,e[1],e[2])
			e[0] = 'done';
	request.messages = filter(lambda x: x[0]!='done',evt);
	#print "Filtered: " + str(request.messages);
	# Now create the image from

	img = Image.new( 'RGB', (cfg_width,cfg_height), "white") # create a new black image
	draw = ImageDraw.Draw(img)
#	pixels = img.load()
	# Collect drawables
	X=[];Y=[];
	if 'x'  in request.session:
		X = request.session['x'];
	if 'y'  in request.session:
		Y = request.session['y'];
	print "Len(X): "+str(len(X))
	for i in range(len(X)-1):
		print "putting at i="+str(i);
		#if (X[i] < 0 or X[i] >= cfg_width or Y[i] < 0 or Y[i] >= cfg_height):
		#	continue;
		draw.line((X[i],Y[i],X[i+1],Y[i+1]),fill=(0,0,0));

	del draw

	with tempfile.TemporaryFile(suffix=".png") as tmpfile:
		img.save(tmpfile, "PNG")
		tmpfile.seek(0)
		data = tmpfile.read();
		s = b64encode(tmpfile.read())
		print s
	return data

def handle_click(request,x,y):
	if 'x' not in request.session:
		request.session['x'] = [x]
	request.session['x'].append(x);

	if 'y' not in request.session:
		request.session['y'] = [y]
	request.session['y'].append(y);

	print("Click handled")
