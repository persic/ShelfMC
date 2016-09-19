import os, sys, ROOT, numpy as n, matplotlib.pyplot as plt

#title=sys.argv[1]
title="Moores bay vs SP"
fileheader=title.replace(" ","_")

#infNames=["PointOneMeterR_parsed.root","FourMeterR_parsed.root","TenMeterR_parsed.root","ThirtyMeterR_parsed.root","FiftyMeterR_parsed.root","OneHundredMeterR_parsed.root","TwoHundredMeterR_parsed.root","FiveHundredMeterR_parsed.root"]
#labels=["0.1m","4.0m","10.0m","30.0m","50.0m","100.0m","200.0m","500.0m"]

#infNames=["FourMeterR_parsed.root","TenMeterR_parsed.root","ThirtyMeterR_parsed.root","FiftyMeterR_parsed.root","OneHundredMeterR_parsed.root","TwoHundredMeterR_parsed.root","FiveHundredMeterR_parsed.root"]
#labels=["4.0m","10.0m","30.0m","50.0m","100.0m","200.0m","500.0m"]

infNames=["VeffForAnalysis/parsedROOTfiles/ShadowParsedLogs.root","VeffForAnalysis/parsedROOTfiles/NoShadowParsedLogs.root","/pub/cpersich/ShelfMC/myshelfmc/ARIANNA_at_SP/parsedRootFiles/ARIANNA_at_SP_Single_WithShadow.root","/pub/cpersich/ShelfMC/myshelfmc/ARIANNA_at_SP/parsedRootFiles/ARIANNA_at_SP_Single_NoShadow.root"]
labels=["Baseline Shadowing","Baseline No Shadowing","At SP Shadowing","At SP No Shadowing"]

#infNames=eval(sys.argv[2])
#labels=eval(sys.argv[3])

inf={}
nt={}
Values={}
VvsEAve={}
VvsENuE={}
VvsENuMu={}
VvsENuTau={}
ErrVsEAve={}
ErrVsENuE={}
ErrVsENuMu={}
ErrVsENuTau={}

for i in range(0,len(labels)):
	inf[i]=ROOT.TFile(infNames[i])
	nt[i]=inf[i].Get("nt")
	Values[i]=[]
	for ent in nt[i]:
		E=ent.energy
		VeffAve=ent.VeffAve
		VeffNuE=ent.VeffNuE
		VeffNuMu=ent.VeffNuMu
		VeffNuTau=ent.VeffNuTau
		WtNuE=ent.wtNuE
		WtNuMu=ent.wtNuMu
		WtNuTau=ent.wtNuTau
		WtAve=(WtNuE+WtNuMu+WtNuTau)/3
		Values[i].append((E,VeffAve,VeffNuE,VeffNuMu,VeffNuTau,WtAve,WtNuE,WtNuMu,WtNuTau))

	Values[i].sort()

	VvsEAve[i]={}
	VvsENuE[i]={}
	VvsENuMu[i]={}
	VvsENuTau[i]={}
	ErrVsEAve[i]={}
	ErrVsENuE[i]={}
	ErrVsENuMu[i]={}
	ErrVsENuTau[i]={}

	for (E,VeffAve,VeffNuE,VeffNuMu,VeffNuTau,WtAve,WtNuE,WtNuMu,WtNuTau) in Values[i]:
		if E in VvsEAve[i]:
			#First do a weighted sum of Veff for each energy
			VvsEAve[i][E]+=VeffAve*WtAve
			VvsENuE[i][E]+=VeffNuE*WtNuE
			VvsENuMu[i][E]+=VeffNuMu*WtNuMu
			VvsENuTau[i][E]+=VeffNuTau*WtNuTau
			#First sum all the weights
			ErrVsEAve[i][E]+=WtAve
			ErrVsENuE[i][E]+=WtNuE
			ErrVsENuMu[i][E]+=WtNuMu
			ErrVsENuTau[i][E]+=WtNuTau
		else:
			#First do a weighted sum of Veff for each energy
			VvsEAve[i][E]=VeffAve*WtAve
			VvsENuE[i][E]=VeffNuE*WtNuE
			VvsENuMu[i][E]=VeffNuMu*WtNuMu
			VvsENuTau[i][E]=VeffNuTau*WtNuTau
			#First sum all the weight
			ErrVsEAve[i][E]=WtAve
			ErrVsENuE[i][E]=WtNuE
			ErrVsENuMu[i][E]=WtNuMu
			ErrVsENuTau[i][E]=WtNuTau

	for E in VvsEAve[i]:
		if ErrVsEAve[i][E] !=0 :
			#divide the weighted sum by the total weight to get the weighted average 
			VvsEAve[i][E]=VvsEAve[i][E]/ErrVsEAve[i][E]
			#StdErr = mean/sqrt(N)
			ErrVsEAve[i][E]=VvsEAve[i][E]/n.sqrt(ErrVsEAve[i][E])
		else:
			VvsEAve[i][E]=0

		if ErrVsENuE[i][E] != 0:
			VvsENuE[i][E]=VvsENuE[i][E]/ErrVsENuE[i][E]
			ErrVsENuE[i][E]=VvsENuE[i][E]/n.sqrt(ErrVsENuE[i][E])
		else:
			VvsENuE[i][E]=0

		if ErrVsENuMu[i][E] != 0:
			VvsENuMu[i][E]=VvsENuMu[i][E]/ErrVsENuMu[i][E]
			ErrVsENuMu[i][E]=VvsENuMu[i][E]/n.sqrt(ErrVsENuMu[i][E])
		else:
			VvsENuMu[i][E]=0

		if ErrVsENuTau[i][E] != 0:
			VvsENuTau[i][E]=VvsENuTau[i][E]/ErrVsENuTau[i][E]
			ErrVsENuTau[i][E]=VvsENuTau[i][E]/n.sqrt(ErrVsENuTau[i][E])
		else:
			VvsENuTau[i][E]=0

plt.figure(1,figsize=(12,9))
plt.suptitle(title)

plt.subplot(2,2,1)
for i in range(0,len(labels)):
	plt.errorbar(VvsEAve[i].keys(),VvsEAve[i].values(),yerr=ErrVsEAve[i].values(),label=labels[i],linestyle='None',marker='o')
plt.xlabel("E")
plt.ylabel("VeffAve (km3Str)")
plt.yscale("log")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)

plt.subplot(2,2,2)
for i in range(0,len(labels)):
	plt.errorbar(VvsENuE[i].keys(),VvsENuE[i].values(),yerr=ErrVsENuE[i].values(),label=labels[i],linestyle='None',marker='o')
plt.xlabel("E")
plt.ylabel("VeffNuE (km3Str)")
plt.yscale("log")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)

plt.subplot(2,2,3)
for i in range(0,len(labels)):
	plt.errorbar(VvsENuMu[i].keys(),VvsENuMu[i].values(),yerr=ErrVsENuMu[i].values(),label=labels[i],linestyle='None',marker='o')
plt.xlabel("E (eV)")
plt.ylabel("VeffNuMu (km3Str)")
plt.yscale("log")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)

plt.subplot(2,2,4)
for i in range(0,len(labels)):
	plt.errorbar(VvsENuTau[i].keys(),VvsENuTau[i].values(),yerr=ErrVsENuTau[i].values(),label=labels[i],linestyle='None',marker='o')
plt.xlabel("E")
plt.ylabel("VeffNuTau (km3Str)")
plt.yscale("log")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)
plt.legend(loc='lower right')
#plt.savefig(fileheader+".png",bbox_inches='tight')
plt.show()
