import os, sys, ROOT, numpy as np, matplotlib.pyplot as plt

#title=sys.argv[1]
title="Shadowing Model Comparison at South Pole"
fileheader=title.replace(" ","_")


infNames=["ParsedLogs/VeffStandardShift_SP_ST2_parsed.root","ParsedLogs/VeffHProp0m_SP_ST2_parsed.root","ParsedLogs/VeffHProp400m_SP_ST2_parsed.root","ParsedLogs/VeffHProp500m_SP_ST2_parsed.root","ParsedLogs/VeffHProp800m_SP_ST2_parsed.root","ParsedLogs/VeffHProp2000m_SP_ST2_parsed.root","ParsedLogs/VeffNoShadowShift_SP_ST2_parsed.root"]
labels=["Standard","H_Prop 0m","H_Prop 400m","H_Prop 500m","H_Prop 800m","H_Prop 2000m","No Shadow"]

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
LoErrVsEAve={}
LoErrVsENuE={}
LoErrVsENuMu={}
LoErrVsENuTau={}
HiErrVsEAve={}
HiErrVsENuE={}
HiErrVsENuMu={}
HiErrVsENuTau={}
RatioVvsEAve={}
RatioVvsENuE={}
RatioVvsENuMu={}
RatioVvsENuTau={}


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
	LoErrVsEAve[i]={}
	LoErrVsENuE[i]={}
	LoErrVsENuMu[i]={}
	LoErrVsENuTau[i]={}
	HiErrVsEAve[i]={}
	HiErrVsENuE[i]={}
	HiErrVsENuMu[i]={}
	HiErrVsENuTau[i]={}
	RatioVvsEAve[i]={}
	RatioVvsENuE[i]={}
	RatioVvsENuMu[i]={}
	RatioVvsENuTau[i]={}

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
			ErrVsEAve[i][E]=VvsEAve[i][E]/np.sqrt(ErrVsEAve[i][E])
		else:
			VvsEAve[i][E]=0

		if ErrVsENuE[i][E] != 0:
			VvsENuE[i][E]=VvsENuE[i][E]/ErrVsENuE[i][E]
			ErrVsENuE[i][E]=VvsENuE[i][E]/np.sqrt(ErrVsENuE[i][E])
		else:
			VvsENuE[i][E]=0

		if ErrVsENuMu[i][E] != 0:
			VvsENuMu[i][E]=VvsENuMu[i][E]/ErrVsENuMu[i][E]
			ErrVsENuMu[i][E]=VvsENuMu[i][E]/np.sqrt(ErrVsENuMu[i][E])
		else:
			VvsENuMu[i][E]=0

		if ErrVsENuTau[i][E] != 0:
			VvsENuTau[i][E]=VvsENuTau[i][E]/ErrVsENuTau[i][E]
			ErrVsENuTau[i][E]=VvsENuTau[i][E]/np.sqrt(ErrVsENuTau[i][E])
		else:
			VvsENuTau[i][E]=0

for i in range(0,len(labels)):
	for E in VvsEAve[i]:
		if VvsEAve[0][E] !=0:
			RatioVvsEAve[i][E]=VvsEAve[i][E]/VvsEAve[0][E]
			HiErrVsEAve[i][E]=(VvsEAve[i][E]+ErrVsEAve[i][E])/(VvsEAve[0][E]-ErrVsEAve[0][E]) - RatioVvsEAve[i][E]
			LoErrVsEAve[i][E]=np.abs((VvsEAve[i][E]-ErrVsEAve[i][E])/(VvsEAve[0][E]+ErrVsEAve[0][E]) - RatioVvsEAve[i][E])
		else:
			RatioVvsEAve[i][E]=0
			HiErrVsEAve[i][E]=0
			LoErrVsEAve[i][E]=0

		if VvsENuE[0][E] !=0:
			RatioVvsENuE[i][E]=VvsENuE[i][E]/VvsENuE[0][E]
			HiErrVsENuE[i][E]=(VvsENuE[i][E]+ErrVsENuE[i][E])/(VvsENuE[0][E]-ErrVsENuE[0][E]) - RatioVvsENuE[i][E]
			LoErrVsENuE[i][E]=np.abs((VvsENuE[i][E]-ErrVsENuE[i][E])/(VvsENuE[0][E]+ErrVsENuE[0][E]) - RatioVvsENuE[i][E])
		else:
			RatioVvsENuE[i][E]=0
			HiErrVsENuE[i][E]=0
			LoErrVsENuE[i][E]=0

		if VvsENuMu[0][E] !=0:
			RatioVvsENuMu[i][E]=VvsENuMu[i][E]/VvsENuMu[0][E]
			HiErrVsENuMu[i][E]=(VvsENuMu[i][E]+ErrVsENuMu[i][E])/(VvsENuMu[0][E]-ErrVsENuMu[0][E]) - RatioVvsENuMu[i][E]
			LoErrVsENuMu[i][E]=np.abs((VvsENuMu[i][E]-ErrVsENuMu[i][E])/(VvsENuMu[0][E]+ErrVsENuMu[0][E]) - RatioVvsENuMu[i][E])
		else:
			RatioVvsENuMu[i][E]=0
			HiErrVsENuMu[i][E]=0
			LoErrVsENuMu[i][E]=0

		if VvsENuTau[0][E] !=0:
			RatioVvsENuTau[i][E]=VvsENuTau[i][E]/VvsENuTau[0][E]
			HiErrVsENuTau[i][E]=(VvsENuTau[i][E]+ErrVsENuTau[i][E])/(VvsENuTau[0][E]-ErrVsENuTau[0][E]) - RatioVvsENuTau[i][E]
			LoErrVsENuTau[i][E]=np.abs((VvsENuTau[i][E]-ErrVsENuTau[i][E])/(VvsENuTau[0][E]+ErrVsENuTau[0][E]) - RatioVvsENuTau[i][E])
		else:
			RatioVvsENuTau[i][E]=0
			HiErrVsENuTau[i][E]=0
			LoErrVsENuTau[i][E]=0

plt.figure(1,figsize=(12,9))
plt.suptitle(title)


plt.subplot(2,2,1)
for i in range(0,len(labels)):
	plt.errorbar(RatioVvsEAve[i].keys(),RatioVvsEAve[i].values(),yerr=[LoErrVsEAve[i].values(),HiErrVsEAve[i].values()],label=labels[i],fmt='o')
plt.xlabel("E")
plt.ylabel("Ratio Ave")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)
plt.ylim(ymin=0.0)

plt.subplot(2,2,2)
for i in range(0,len(labels)):
	plt.errorbar(RatioVvsENuE[i].keys(),RatioVvsENuE[i].values(),yerr=[LoErrVsENuE[i].values(),HiErrVsENuE[i].values()],label=labels[i],fmt='o')
plt.xlabel("E")
plt.ylabel("Ratio NuE")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)
plt.ylim(ymin=0.0)

plt.subplot(2,2,3)
for i in range(0,len(labels)):
	plt.errorbar(RatioVvsENuMu[i].keys(),RatioVvsENuMu[i].values(),yerr=[LoErrVsENuMu[i].values(),HiErrVsENuMu[i].values()],label=labels[i],fmt='o')
plt.xlabel("E (eV)")
plt.ylabel("Ratio NuMu")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)
plt.ylim(ymin=0.0)

plt.subplot(2,2,4)
for i in range(0,len(labels)):
	plt.errorbar(RatioVvsENuTau[i].keys(),RatioVvsENuTau[i].values(),yerr=[LoErrVsENuTau[i].values(),HiErrVsENuTau[i].values()],label=labels[i],fmt='o')
plt.xlabel("E")
plt.ylabel("Ratio NuTau")
plt.xscale("log")
plt.gca().xaxis.grid(True)
plt.gca().yaxis.grid(True)
plt.ylim(ymin=0.0)
plt.legend()
#plt.savefig(fileheader+".png",bbox_inches='tight')
plt.show()
