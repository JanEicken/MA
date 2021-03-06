import trajcomp;
import sys, getopt, os;
import numpy
import math

def curvature(b) :
    #a = numpy.array([ [  0.  ,   0.  ],[  0.3 ,   0.  ],[  1.25,  -0.1 ],
    #          [  2.1 ,  -0.9 ],[  2.85,  -2.3 ],[  3.8 ,  -3.95],
    #          [  5.  ,  -5.75],[  6.4 ,  -7.8 ],[  8.05,  -9.9 ],
    #          [  9.9 , -11.6 ],[ 12.05, -12.85],[ 14.25, -13.7 ],
    #          [ 16.5 , -13.8 ],[ 19.25, -13.35],[ 21.3 , -12.2 ],
    #          [ 22.8 , -10.5 ],[ 23.55,  -8.15],[ 22.95,  -6.1 ],
    #          [ 21.35,  -3.95],[ 19.1 ,  -1.9 ]])
    c = []
    for i in b:
        c.append([i[0]*(numpy.pi/180), i[1]*(numpy.pi/180)])
    print c
    #print "Do another one"
    a = numpy.array(c)
    dx_dt = numpy.gradient(a[:, 0])
    dy_dt = numpy.gradient(a[:, 1])
    ds_dt = numpy.sqrt(dx_dt * dx_dt + dy_dt * dy_dt)
    d2s_dt2 = numpy.gradient(ds_dt)
    d2x_dt2 = numpy.gradient(dx_dt)
    d2y_dt2 = numpy.gradient(dy_dt)
    curvature = numpy.abs(d2x_dt2 * dy_dt - dx_dt * d2y_dt2) / (dx_dt * dx_dt + dy_dt * dy_dt)**1.5
    curv = curvature 
    print "CURVS"
    print curv
    return curv

def curvatur(b):
    2*((x2-x1)*(y3-y2)-(y2-y1)*(x3-x2)) / sqrt(
        ((x2-x1)^2+(y2-y1)^2)*((x3-x2)^2+(y3-y2)^2)*((x1-x3)^2+(y1-y3)^2) );

def calcCurvature(coords_file):
    h = 0.000001
    res = []
    with open(coords_file) as f:
        lines = f.read().splitlines()
        curvs = []

        for l in range(0, len(lines)):
            #print "Line::"
            print "line"
            print lines[l]
            traj = lines[l].split(' ')

            res.append([float(traj[0]), float(traj[1])])
            continue
            if l == 0 or l == len(lines)-1:
                curv = 0.0
                continue
            prev = lines[l-1].split(' ');
            nex = lines[l+1].split(' ');


            xCoefficientArray = [float(traj[0]) - float(prev[0]), float(traj[1]) - float(prev[1])]
            yCoefficientArray = [float(nex[0]) - float(prev[0]), float(nex[1]) - float(prev[1])]

            coefficientArray = numpy.array([xCoefficientArray,yCoefficientArray])
            constantArray = numpy.array([(pow(float(traj[0]),2) + pow(float(traj[1]),2) - pow(float(prev[0]),2) - pow(float(prev[1]),2))/2, (pow(float(nex[0]),2) + pow(float(nex[1]),2) - pow(float(prev[0]),2) - pow(float(prev[1]),2))/2])
            try:
                center = numpy.linalg.solve(coefficientArray, constantArray)

                #print "center: " + str(center)
                radius_prev = math.sqrt(pow(float(prev[0])-center[0],2) + pow(float(prev[1])-center[1],2))
                radius_cur = math.sqrt(pow(float(traj[0])-center[0],2) + pow(float(traj[1])-center[1],2))
                radius_nex = math.sqrt(pow(float(nex[0])-center[0],2) + pow(float(nex[1])-center[1],2))

                #print "Radius: " + str(radius)
                curv = 1/radius_prev
                print "last curv: " + str((1/radius_prev)% 360.0)+ " cur curv " + str((1/radius_cur)% 360.0)  + " nex curv" + str((1/radius_nex)% 360.0)
                #print "curv: " + str(curv)
                #print math.sqrt(pow(p_prev[0]-center[0],2) + pow(p_prev[1]-center[1],2))
                #print math.sqrt(pow(float(traj[0])-center[0],2) + pow(float(traj[1])-center[1],2))
                #print math.sqrt(pow(float(nex[0])-center[0],2) + pow(p_next[1]-center[1],2))

            except:
                print "Straight line"
                curv = 0.0
            curvs.append(curv)

            continue


            traj_last = lines[l-1].split(' ')
            traj_next = lines[l+1].split(' ')
            #print "Cur point " + ', '.join(traj)
            #print "Last point " + ', '.join(traj_last)
            #print "Next point " + ', '.join(traj_next)
            d_prev = math.sqrt(pow(float(traj[0])-float(traj_last[0]),2) + pow(float(traj[1])-float(traj_last[1]),2))
            d_next = math.sqrt(pow(float(traj_next[0])-float(traj[0]),2) + pow(float(traj[1])-float(traj_next[1]),2))
            if(d_prev == 0) or d_next == 0:
                curv = 0.0
                continue
            d_prev = d_prev*1000
            d_next = d_next*1000
            #print "distance prev: " + str(d_prev)
            #print "distance next: " + str(d_next)
            p_prev = []
            p_next = []

            if(float(traj[0]) < float(traj_last[0])):
                p_prev.append(float(traj[0]) + float(traj[0])*(h/d_prev));
            else:
                p_prev.append(float(traj[0]) - float(traj[0])*(h/d_prev));

            if(float(traj[1]) < float(traj_last[1])):
                p_prev.append(float(traj[1]) + float(traj[1])*(h/d_prev));
            else:
                p_prev.append(float(traj[1]) - float(traj[1])*(h/d_prev));

            if(float(traj_next[0]) < float(traj[0])):
                p_next.append(float(traj[0]) + float(traj[0])*(h/d_next));
            else:
                p_next.append(float(traj[0]) - float(traj[0])*(h/d_next));

            if(float(traj_next[1]) < float(traj[1])):
                p_next.append(float(traj[1]) + float(traj[1])*(h/d_next));
            else:
                p_next.append(float(traj[1]) - float(traj[1])*(h/d_next));

            #print "last calculated point: " + str(p_prev[0]) + " " + str(p_prev[1])
            #print "next calculated point: " + str(p_next[0]) + " " + str(p_next[1])

            curv = 0.0;
            #print p_prev[0][2:]
            xCoefficientArray = [float(traj[0]) - p_prev[0], float(traj[1]) - p_prev[1]]
            yCoefficientArray = [p_next[0] - p_prev[0], p_next[1] - p_prev[1]]

            coefficientArray = numpy.array([xCoefficientArray,yCoefficientArray])
            constantArray = numpy.array([(pow(float(traj[0]),2) + pow(float(traj[1]),2) - pow(p_prev[0],2) - pow(p_prev[1],2))/2, (pow(p_next[0],2) + pow(p_next[1],2) - pow(p_prev[0],2) - pow(p_prev[1],2))/2])
            try:
                center = numpy.linalg.solve(coefficientArray, constantArray)

                #print "center: " + str(center)
                radius = math.sqrt(pow(p_prev[0]-center[0],2) + pow(p_prev[1]-center[1],2))
                #print "Radius: " + str(radius)
                curv = 1/radius
                #print "curv: " + str(curv)
                #print math.sqrt(pow(p_prev[0]-center[0],2) + pow(p_prev[1]-center[1],2))
                #print math.sqrt(pow(float(traj[0])-center[0],2) + pow(float(traj[1])-center[1],2))
                #print math.sqrt(pow(p_next[0]-center[0],2) + pow(p_next[1]-center[1],2))

            except:
                print "Straight line"
                curv = 0.0
            curvs.append(curv)
    print "to calculate: "
    print res
    print "---------------------------------------------------------"
    print curvs
    print "---------------------------------------------------------"
    curvature(res)
    return curvs;



def main(argv):

    inputfolder = ''
    outputfolder = ''
    beta = 0
    try:
        opts, args = getopt.getopt(argv,"hi:o:b")
    except getopt.GetoptError:
        print 'persistanceAlgo.py -i <inputfolder> -o <outputfolder> -b <beta>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'persistanceAlgo.py -i <inputfolder> -o <outputfolder> -b <beta>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfolder = arg
        elif opt in ("-o", "--ofile"):
            outputfolder = arg
        elif opt in("-b"):
            beta = arg

    calcCurvature("/home/jan/Desktop/Masterarbeit/trajcomputing_code/v_2/libtrajcomp-src/python/TestTraj.csv")
    return;
    if not os.path.exists(inputfolder):
        print "inputfolder does not exist"
        return;

    print 'Input folder is "', inputfolder
    print 'Output folder is "', outputfolder
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    for i in os.listdir(inputfolder):
        if(i.endswith(".csv")):
            newfolder = i.split(".")
            print "/home/jan/Desktop/Masterarbeit/trajcomputing_code/v_2/libtrajcomp-src/doc/md/data/processed/GoogleNow/testing/"+i
            calcCurvature("/home/jan/Desktop/Masterarbeit/trajcomputing_code/v_2/libtrajcomp-src/doc/md/data/processed/GoogleNow/testing/"+i)

            continue
        else:
            continue

if __name__ == "__main__":
    main(sys.argv[1:])
