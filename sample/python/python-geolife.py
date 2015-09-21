import trajcomp;
from operator import itemgetter 
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

print ("Loading first set of GeoLife Dataset");
A=trajcomp.geolife(15,"/home/martin/datasets/geolife/")
handles = range(A,trajcomp.size());
print "Found " + str(trajcomp.size())
xlimit = [39.689602,40.122410]
ylimit = [116.105789,116.670021]


# list of lists is suboptimal, but does not depend on numpy
for i in handles:
  A=trajcomp.get(i);
  trajcomp.summary(i);
  X = map(itemgetter(0), A) ;
  Y = map(itemgetter(1),A);
  #plt.clf()
  plt.plot(X, Y, '-')
 # Grossraum
#  plt.xlim([39.75,40.05]);
#  plt.ylim([116.116,116.667]);
  plt.xlim(xlimit);
  plt.ylim(ylimit)
  #plt.show()
  plt.savefig(("frame-%04d.png" % i), bbox_inches='tight')
