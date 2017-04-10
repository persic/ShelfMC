import os,sys
import subprocess

def runVeffShelfMC(workDir, outTag):
    args = [ "/pub/cpersich/ShelfMC/myshelfmc/shelfmc_stripped.exe",
             workDir, outTag ]
    print " ".join(args)
    subprocess.call(args)
    

if __name__=="__main__":
    if (len(sys.argv)!=2):
        print "Usage python runVeffShelfMC.py [job file]."
        print ""
        print "Where the job file contains one line with:"
        print "   [working directory]  [output file tag] "
        sys.exit(1)
        
    with open(sys.argv[1]) as inf:
        x = inf.readline().strip().split()
        if (len(x)!=2):
            print "File [{0}] contains an invalid number of parameters [{1}]."\
                .format(sys.argv[1], len(x))
            sys.exit(1)
        else:
            runVeffShelfMC(x[0], x[1])

