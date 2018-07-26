import os, sys, ROOT, numpy as np
from plotVeffUtils import *

infn = sys.argv[1]
outfn = sys.argv[2]

if len(sys.argv) != 3:
    print('USAGE: python makeVeffCsv [input filename] [output filename]')
    sys.exit()

DArray=CalcVeffMutiLog(infn)
SaveCSV(DArray,outfn)
