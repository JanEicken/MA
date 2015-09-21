import trajcomp;
from operator import itemgetter 
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


A=trajcomp.load("../../data/prague1.dat");
print "Reference is " + str(A)

# list of lists is suboptimal, but does not depend on numpy
for i in range(1,150):
  epsilon = i/10.0;
  epsilon = epsilon * epsilon;
  print epsilon;
  B = trajcomp.douglas_peucker(A,epsilon);
  trajcomp.summary(B);
  C =trajcomp.get(B)
  
  X = map(itemgetter(0), C) ;
  Y = map(itemgetter(1),C);
  plt.clf();
  fig, ax1 = plt.subplots()
  ax1.plot(X, Y, '-o')
  ax1.set_xlabel('X')
  ax1.set_ylabel('Y')
  
  ax2 = ax1.twinx()
  plt.ylim(ymax=(15*15),ymin=0);
  space=abs(max(X)-min(X))/20;
  
  ax2.bar(max(X)+2*space,epsilon,width=space);
  
  #plt.show()
  plt.savefig(("frame-%04d.png" % i), bbox_inches='tight')
