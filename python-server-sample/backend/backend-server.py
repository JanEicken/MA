import SocketServer
import trajcomp;
import json;
import sys;
import os, shutil

cfg_numdir = 2;
cfg_location = "/home/jan/Desktop/Masterarbeit/trajcomputing_code/v_2/libtrajcomp-src/data/Geolife";
gN_location = "/home/jan/Desktop/Masterarbeit/trajcomputing_code/v_2/libtrajcomp-src/doc/md/data/GoogleNow";
gN_numdir = 1;
elki_location = "/home/jan/elki.jar";
resample_location = "/home/jan/Desktop/Masterarbeit/trajcomputing_code/v_2/libtrajcomp-src/doc/md/resample/"



def top_k_jaccard_id(theid):
	s = trajcomp.jaccard_id(theid);
	print min(s)
	print max(s)

	# Retrieve top k
	k=10
	# First sort and keep indizes
	i = [i[0] for i in sorted(enumerate(s), key=lambda x:x[1])]

	i = i[:k]  #[-k,:]
	for idx in i:
		print "Taking "+str(s[idx]);

	return i;

def top_k_jaccard(sl):
	s = trajcomp.jaccard_stringlist(sl);
	print min(s)
	print max(s)

	# Retrieve top k
	k=10
	# First sort and keep indizes
	i = [i[0] for i in sorted(enumerate(s), key=lambda x:x[1])]

	print i[0]
	i = i[:k]  #[-k,:]
	print i[0]
	for idx in i:
		print "Taking "+str(s[idx]);

	return i

def top_k_intersect(sl):
	s = trajcomp.intersect_stringlist(sl);
	print min(s)
	print max(s)

	tau = min(s) + 0.75 *(max(s)-min(s));

	# Retrieve top k
	k=10
	# First sort and keep indizes
	i = [i[0] for i in sorted(enumerate(s), key=lambda x:x[1])]
#	i = i[-k:]
	r=[];
	for idx in i:
		if(s[idx] > tau):
		  r.append(idx);
		  print "Left "+str(s);
	return r


def top_k_subset_correlation(sl):
	s = trajcomp.relative_subset_contradicts_stringlist(sl);
#	print s
	# subset correlation: take only zero, but double == zero is fabs
	i=[i for i, j in enumerate(s) if j <=10E-6 ]
	# Retrieve up to
#	print i
	print "Found " +str(len(i)) + " fully correlating. Limit to k=10"
	k=10
	i = i[:k]  #[-k,:]
	sys.stdout.flush()
	return i

def dbscan_segmentation(d, m, e):
	print  "dbscan backend server";
	print d;
	print m;
	print e;
	print elki_location;
	r = trajcomp.dbscan(float(e),int(m), "Test", d, elki_location);
	return r;

def douglas_peucker(d, e):
	r = trajcomp.douglas_peucker_online(d,float(e));
	return r;

def resample_traj(d, h):
	r = trajcomp.resample(d, h, resample_location);
	print "got back:"
	print r
	#Delete created files
	#for the_file in os.listdir(resample_location):
	#    file_path = os.path.join(resample_location, the_file)
	#    try:
	#        if os.path.isfile(file_path):
	#            os.unlink(file_path)
	#        #elif os.path.isdir(file_path): shutil.rmtree(file_path)
	#    except Exception, e:
	#        print e
	return r;

def getID(i):
	p = "../../doc/md/data/processed/GoogleNow/Jan/" + str(i) + ".csv";
	pt = "../../doc/md/data/processed/GoogleNow/Jan/times/" + str(i) + ".csv";
	r = trajcomp.get_g_id(p, pt)
	print "got back: "

	return r;


class MyTCPHandler(SocketServer.StreamRequestHandler):
    """
    The RequestHandler class for our server.

    It is instantiated once per connection to the server, and must
    override the handle() method to implement communication to the
    client.
    """
    def handle(self):
        # self.rfile is a file-like object created by the handler;
        # we can now use e.g. readline() instead of raw recv() calls
        self.data="";
        while(self.data != "QUIT"):
			self.data = self.rfile.readline().strip()
			C=self.data.split();
			if C[0] == "GET":
				T = trajcomp.get(int(C[1]))
				self.wfile.write(json.dumps(T)+"\n")
			if C[0] == "SUBSETCORR":
				answer = top_k_subset_correlation(C[1:]);
				self.wfile.write(json.dumps(answer)+"\n");
			if C[0] == "JACCARD":
				answer = top_k_jaccard(C[1:]);
				self.wfile.write(json.dumps(answer)+"\n");
			if C[0] == "JACCARDID":
				answer = top_k_jaccard_id(int(C[1]));
				self.wfile.write(json.dumps(answer)+"\n");
			if C[0] == "INTERSECT":
				#print C[1:];
				answer = top_k_intersect(C[1:]);
				self.wfile.write(json.dumps(answer)+"\n");
			if C[0] == "DBSCAN":
				print C[1:];
				answer = dbscan_segmentation(C[1], C[2], C[3]);
				self.wfile.write(json.dumps(answer)+"\n");
			if C[0] == "DOUGLAS":
				print "Backendserver";
				print C[1:];
				T = trajcomp.douglas_peucker_online(C[1],float(C[2]));
				print "Result DP";
				print T;
				self.wfile.write(json.dumps(T)+"\n");
			if C[0] == "RESAMPLE":
				print "Resample backend";
				print C[1:];
				answer = resample_traj(C[1], int(C[2]));
				self.wfile.write(json.dumps(answer)+"\n");
			if C[0] == "THRESHOLD":
				print "Threshold backend";
				print C[1:];
				T = trajcomp.threshold(float(C[1]), float(C[2]), "../TestTraj.csv", "../TestTrajTimes.csv");
				self.wfile.write(json.dumps(T)+"\n");
			if C[0] == "PERSISTENCE":
				print "Peristence backend";
				print C[1:];
				T = trajcomp.persistence(float(C[1]), "../TestTraj.csv", "../TestTrajTimes.csv");
				self.wfile.write(json.dumps(T)+"\n");
			if C[0] == "GETGID":
				print "GETID backend";
				print C[1:];
				T = getID(C[1]);
				self.wfile.write(json.dumps(T)+"\n");

			if self.data=="HELO":
				self.wfile.write("HELO too\n");
			if self.data=="KILL":
				self.server.shutdown();
				self.data="QUIT";


class ThreadingServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer): pass

if __name__ == "__main__":
	HOST, PORT = "localhost", 9999

    # Create the server, binding to localhost on port 9999
	server = ThreadingServer((HOST, PORT), MyTCPHandler)

	print ("Loading first set of GeoLife Dataset");
	print cfg_location;
	#A=trajcomp.geolife(cfg_numdir, cfg_location)
	# handles = range(A,trajcomp.size());
	#print "Found " + str(trajcomp.size())
	# Create the index
	#trajcomp.create_bloomindex();

	server.serve_forever()
