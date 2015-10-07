import socket
import json;

# GET A TRAJECTORY FROM PERSISTENT BACKEND MODULE

def backend_get_g_id(i):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	fss += i + " ";
	fss = fss[:-1]
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("GETGID "+fss+" \n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	print "backendclient back"
	return json.loads(s);

def backend_threshold_data(p):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	fss += p + " ";
	fss = fss[:-1]
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("PERSISTENCE "+fss+" \n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

def backend_threshold_data(v, o):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	fss += v + " " + o + " ";
	fss = fss[:-1]
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("THRESHOLD "+fss+" \n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

def backend_resample_data(d, h):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	fss += d + " " + h + " ";
	fss = fss[:-1]
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("RESAMPLE "+fss+" \n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

def backend_DP_segmentation(d, eps):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	fss += d + " " + eps + " ";
	fss = fss[:-1]
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("DOUGLAS "+fss+"\n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

def backend_DBSCAN_segmentation(d, minPts, eps):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	fss += d + " " + minPts + " " +  eps + " ";
	fss = fss[:-1]

	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("DBSCAN "+fss+"\n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

def backend_get(i):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("GET "+str(i)+"\n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

# GET SIMILAR TO ITEM

def backend_getjaccard(idx):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("JACCARDID "+idx+"\n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);


# GET A LIST OF INDICES BASED ON FILTERSET STRINGS

def backend_jaccard(filterset):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	for x in filterset:
		fss += x + " ";

	fss = fss[:-1];
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("JACCARD "+fss+"\n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

# GET A LIST OF INDICES BASED ON FILTERSET STRINGS

def backend_intersect(filterset):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	for x in filterset:
		fss += x + " ";

	fss = fss[:-1];
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("INTERSECT "+fss+"\n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);

def backend_subsetcorr(filterset):
	HOST, PORT = "localhost", 9999
	# Create a socket (SOCK_STREAM means a TCP socket)
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	fss = "";
	for x in filterset:
		fss += x + " ";

	fss = fss[:-1];
	try:
		sock.connect((HOST, PORT))
		f = sock.makefile()
		sock.sendall("SUBSETCORR "+fss+"\n");
		#received = sock.recv(1024)
		s = f.readline();
		sock.sendall("QUIT\n");
	finally:
		sock.close()
	return json.loads(s);


#print backend_intersect(["abc","def"]);
