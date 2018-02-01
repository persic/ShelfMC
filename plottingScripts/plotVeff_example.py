import os, sys, ROOT, numpy as n, matplotlib.pyplot as plt
from plotVeffUtils import *

title="Moores bay vs SP"
fileheader=title.replace(" ","_")
infNames=["ParsedLogs_Example/VeffCombo_3Sig_MB_HProp500_ST6_parsed.root","ParsedLogs_Example/VeffCombo_3Sig_SP_1400m_HProp500_ST6_parsed.root"]
labels=["Moore's Bay","South Pole"]

NumFiles=len(labels)

plt.figure(1,figsize=(12,9))
plt.suptitle(title)
ax1=plt.subplot(2,2,1)
ax2=plt.subplot(2,2,2)
ax3=plt.subplot(2,2,3)
ax4=plt.subplot(2,2,4)
#DArrays=[]

for i in range(0,NumFiles):
	DArray=CalcVeffMutiLog(infNames[i])
	#DArrays.append(DArray)
	E, VeffAve, VeffNuE, VeffNuMu, VeffNuTau, VeffErrAve, VeffErrNuE, VeffErrNuMu, VeffErrNuTau = GetAllColumns(DArray)

	ax1.errorbar(E,VeffAve,yerr=VeffErrAve,ls='None',marker='o',label=labels[i])
	ax1.set_xlabel("E (eV)")
	ax1.set_ylabel("VeffAve (km^3Str)")
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.grid()

	ax2.errorbar(E,VeffNuE,yerr=VeffErrNuE,ls='None',marker='o',label=labels[i])
	ax2.set_xlabel("E (eV)")
	ax2.set_ylabel("VeffNuE (km^3Str)")
	ax2.set_xscale("log")
	ax2.set_yscale("log")
	ax2.grid()

	ax3.errorbar(E,VeffNuMu,yerr=VeffErrNuMu,ls='None',marker='o',label=labels[i])
	ax3.set_xlabel("E (eV)")
	ax3.set_ylabel("VeffNuMu (km^3Str)")
	ax3.set_xscale("log")
	ax3.set_yscale("log")
	ax3.grid()

	ax4.errorbar(E,VeffNuTau,yerr=VeffErrNuTau,ls='None',marker='o',label=labels[i])
	ax4.set_xlabel("E (eV)")
	ax4.set_ylabel("VeffNuTau (km^3Str)")
	ax4.set_xscale("log")
	ax4.set_yscale("log")
	ax4.grid()

ax4.legend(loc=4)
plt.show()
