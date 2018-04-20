import sys, os
import ROOT
import root_numpy
import numpy as np
from EventListTools import *

infn = sys.argv[1]
EventListinfn = sys.argv[2]
outfn = sys.argv[3]

Fin = ROOT.TFile(infn)
PAM = Fin.Get("PAM")

EvtAndWeight =  root_numpy.tree2array(PAM,['ievt','weight'])

header = 'evid weight'
typeFormat = ['%i','%f']

Events = ReadEventListARASimToListObject(EventListinfn)

weights = np.zeros(len(Events))
triggered = np.zeros(len(Events))

OutArray = ARASimListArray(Events)

for i, w in EvtAndWeight:
    triggered[i] = 1
    weights[i] = w

OutArray = np.column_stack((OutArray,weights,triggered))

header = 'evid nuflavorint nu_nubar pnu currentint posnu_r posnu_theta posnu_phi nnu_theta nnu_phi elast_y weight globalTrig'
typeFormat = ['%i','%i','%i','%f','%i','%f','%f','%f','%f','%f','%f','%f','%i']


np.savetxt(outfn,OutArray,header=header,fmt=typeFormat,delimiter=' ')
