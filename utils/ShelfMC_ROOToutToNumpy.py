import sys, os
import ROOT
import root_numpy
import numpy as np
from EventListTools import *

infn = sys.argv[1]
outfn = sys.argv[2]

Fin = ROOT.TFile(infn)
PAM = Fin.Get("PAM")

ShelfMCToPDGParticleConvert = {
    1:12,
    2:14,
    3:16
    }

ParamList = {'evtNum':'ievt',
            'log10E':'Energy',
            'flavor':'flavor',
            'current':'mycurrent',
            'xyz':'Posi_Int',
            'inelasticity':'y',
            'eshower_em':'eshower_em',
            'eshower_had':'eshower_had',
            'arrivalTheta':'theta_signal_atAT',
            'arrivalTheta_Ant':'theta_my_signalAT',
            'arrivalPhi':'phi_nposnu2ST',
            'nnu_theta':'theta_nu',
            'nnu_phi':'phi_nu',
            'viewangle':'viewangle_triggered',
            #'viewangle_Ant':'viewangle_triggered_LPA',
            'coneAngle':'changle',
            'n_pol':'Fresnel_Pol_D',
            'freq':'freq_my',
            'distance':'dis',
            'attenfactor':'attenfactor',
            #'EField':'vmmhz_my',
            'vMaxTrig':'volte_allAT_mymax',
            'vTrig_Ant':'volt_LPA',
            'weight':'weight',
            'nStnTrigDirect':'sum_triggeredST',
            'nStnTrigRefl':'sum_triggeredST_mirror'
            }

ParameterArray =  root_numpy.tree2array(PAM,ParamList.values())
ParameterArray.dtype.names = tuple(ParamList.keys())

np.save(outfn,ParameterArray)
