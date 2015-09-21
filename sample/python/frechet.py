import trajcomp;
import sys;
from operator import itemgetter 
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import time;

def elapsed(var):
	t = time.clock();
	r = eval(var)
	return (r,time.clock()-t)

epsilon = float(sys.argv[1]);

print ("Loading first set of GeoLife Dataset");
A=trajcomp.geolife(1,"/home/martin/datasets/geolife/")
handles = range(A,trajcomp.size());
print "Found " + str(trajcomp.size())

# list of lists is suboptimal, but does not depend on numpy
for i in handles:
	for j in handles:
		(v_real,t_real) = elapsed("trajcomp.frechet(%d,%d,-1)" % (i,j))
		A = trajcomp.douglas_peucker(i,epsilon)
		B = trajcomp.douglas_peucker(j,epsilon)
#		trajcomp.summary(i);
#		trajcomp.summary(j);
#		trajcomp.summary(A);
#		trajcomp.summary(B);
		(v_sim,t_sim) = elapsed("trajcomp.frechet(%d,%d,-1)" % (A,B))
		print str(epsilon) + " " + str(i) +" "+ str(j)+" " + str (t_real) + " " + str(t_sim) + " "  +str(v_real)+" " + 	str(v_sim) +" " + str(v_real / v_sim)
        	
		
		
		
