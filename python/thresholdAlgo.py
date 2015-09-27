import trajcomp;
import sys, getopt, os;

elki_location = "/home/jan/elki.jar";
def main(argv):
    inputfolder = ''
    outputfolder = ''
    orientationThresh = 0
    velocityThresh = 0
    try:
        opts, args = getopt.getopt(argv,"hi:o:t:v:")
    except getopt.GetoptError:
        print 'test.py -i <inputfolder> -o <outputfolder> -t <orientation Threshold> -v <velocity Threshold>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'densityAlgo.py -i <inputfolder> -o <outputfolder> -t <orientation Threshold> -v <velocity Threshold>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfolder = arg
        elif opt in ("-o", "--ofile"):
            outputfolder = arg
        elif opt in("-t"):
            eps = arg
        elif opt in("-v"):
            minPts = arg

    if not os.path.exists(inputfolder):
        print "inputfolder does not exist"
        return;

    print 'Input folder is "', inputfolder
    print 'Output folder is "', outputfolder
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    for i in os.listdir(inputfolder):
        if(i.endswith(".csv")):
            print i
            newfolder = i.split(".")
            R = trajcomp.threshold(float(velocityThresh), float(orientationThresh),inputfolder+i,
            inputfolder+"times/"+i)
            # Do plots on result
            print R
            continue

        else:
            continue


if __name__ == "__main__":
    main(sys.argv[1:])
