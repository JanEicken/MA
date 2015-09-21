import trajcomp;
from operator import itemgetter 
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def printfig(name,traj):
	C=trajcomp.get(traj);
	X = map(itemgetter(0), C) ;
	Y = map(itemgetter(1), C);
	plt.clf();
	fig, ax1 = plt.subplots()
	ax1.plot(X, Y, '-o')
	ax1.set_xlabel('X')
	ax1.set_ylabel('Y')
	#plt.show()
	plt.savefig(name)


A=trajcomp.load("../../data/ngon10.dat");

B=trajcomp.uniform_select(A,2,True);
C=trajcomp.uniform_select(A,3,True);

printfig("uniform_select_in.png",A);
printfig("uniform_select_2.png",B);
printfig("uniform_select_3.png",C);


