import os, sys, ROOT, numpy as np, matplotlib.pyplot as plt

def DictToArray(Dict):
	outArray=np.array(Dict.items())
	outArray=np.sort(outArray,axis=0)
	return outArray

def CalcVeffMutiLog(infn,UseWeights=True,Tname="nt"):
	inf=ROOT.TFile(infn)
	nt=inf.Get(Tname)
	Values=[]
	for ent in nt:
		E=ent.energy
		VeffAve=ent.VeffAve
		VeffNuE=ent.VeffNuE
		VeffNuMu=ent.VeffNuMu
		VeffNuTau=ent.VeffNuTau
		WtNuE=ent.wtNuE
		WtNuMu=ent.wtNuMu
		WtNuTau=ent.wtNuTau
		WtAve=(WtNuE+WtNuMu+WtNuTau) #This should be the total weighted number of neutrinos
		if UseWeights==False:
			nNuEtrg = ent.nNuEtrg
			nNuMutrg = ent.nNuMutrg
			nNuTautrg = ent.nNuTautrg
			nAvetrg = (nNuEtrg + nNuMutrg + nNuTautrg)/3.0
			if WtNuE != 0:
				VeffNuE = VeffNuE * nNuEtrg / WtNuE
			else:
				VeffNuE = 0
			if WtNuMu != 0:
				VeffNuMu = VeffNuMu * nNuMutrg / WtNuMu
			else:
				VeffNuMu = 0
			if WtNuTau != 0:
				VeffNuTau = VeffNuTau * nNuTautrg / WtNuTau
			else:
				VeffNuTau = 0
			if WtNuE != 0:
				VeffAve = VeffAve * nAvetrg / WtNuE
			else:
				VeffAve = 0
			WtNuE = nNuEtrg
			WtNuMu = nNuMutrg
			WtNuTau = nNuTautrg
			WtAve = nAvetrg
		Values.append((E,VeffAve,VeffNuE,VeffNuMu,VeffNuTau,WtAve,WtNuE,WtNuMu,WtNuTau))
	VvsEAve={}
	VvsENuE={}
	VvsENuMu={}
	VvsENuTau={}
	ErrVsEAve={}
	ErrVsENuE={}
	ErrVsENuMu={}
	ErrVsENuTau={}
	for (E,VeffAve,VeffNuE,VeffNuMu,VeffNuTau,WtAve,WtNuE,WtNuMu,WtNuTau) in Values:
		if E in VvsEAve:
			#First do a weighted sum of Veff for each energy
			VvsEAve[E]+=VeffAve*WtAve
			VvsENuE[E]+=VeffNuE*WtNuE
			VvsENuMu[E]+=VeffNuMu*WtNuMu
			VvsENuTau[E]+=VeffNuTau*WtNuTau
			#First sum all the weights
			ErrVsEAve[E]+=WtAve
			ErrVsENuE[E]+=WtNuE
			ErrVsENuMu[E]+=WtNuMu
			ErrVsENuTau[E]+=WtNuTau
		else:
			#First do a weighted sum of Veff for each energy
			VvsEAve[E]=VeffAve*WtAve
			VvsENuE[E]=VeffNuE*WtNuE
			VvsENuMu[E]=VeffNuMu*WtNuMu
			VvsENuTau[E]=VeffNuTau*WtNuTau
			#First sum all the weight
			ErrVsEAve[E]=WtAve
			ErrVsENuE[E]=WtNuE
			ErrVsENuMu[E]=WtNuMu
			ErrVsENuTau[E]=WtNuTau

	for E in VvsEAve:
		if ErrVsEAve[E] !=0 :
			#divide the weighted sum by the total weight to get the weighted average
			VvsEAve[E]=VvsEAve[E]/ErrVsEAve[E]
			#StdErr = mean/sqrt(N)
			ErrVsEAve[E]=VvsEAve[E]/np.sqrt(ErrVsEAve[E])
		else:
			VvsEAve[E]=0

		if ErrVsENuE[E] != 0:
			VvsENuE[E]=VvsENuE[E]/ErrVsENuE[E]
			ErrVsENuE[E]=VvsENuE[E]/np.sqrt(ErrVsENuE[E])
		else:
			VvsENuE[E]=0

		if ErrVsENuMu[E] != 0:
			VvsENuMu[E]=VvsENuMu[E]/ErrVsENuMu[E]
			ErrVsENuMu[E]=VvsENuMu[E]/np.sqrt(ErrVsENuMu[E])
		else:
			VvsENuMu[E]=0

		if ErrVsENuTau[E] != 0:
			VvsENuTau[E]=VvsENuTau[E]/ErrVsENuTau[E]
			ErrVsENuTau[E]=VvsENuTau[E]/np.sqrt(ErrVsENuTau[E])
		else:
			VvsENuTau[E]=0

	VvsEAve=DictToArray(VvsEAve)
	VvsENuE=DictToArray(VvsENuE)
	VvsENuMu=DictToArray(VvsENuMu)
	VvsENuTau=DictToArray(VvsENuTau)
	ErrVsEAve=DictToArray(ErrVsEAve)
	ErrVsENuE=DictToArray(ErrVsENuE)
	ErrVsENuMu=DictToArray(ErrVsENuMu)
	ErrVsENuTau=DictToArray(ErrVsENuTau)

	OutArray=np.column_stack((VvsEAve[:,0],VvsEAve[:,1],VvsENuE[:,1],VvsENuMu[:,1],VvsENuTau[:,1],ErrVsEAve[:,1],ErrVsENuE[:,1],ErrVsENuMu[:,1],ErrVsENuTau[:,1]))
	#names = ['E', 'VeffAve', 'VeffNuE', 'VeffNuMu', 'VeffNuTau', 'VeffErrAve', 'VeffErrNuE', 'VeffErrNuMu', 'VeffErrNuTau']
	return OutArray

def GetColumn(DataArray,name):
	names = ['E', 'VeffAve', 'VeffNuE', 'VeffNuMu', 'VeffNuTau', 'VeffErrAve', 'VeffErrNuE', 'VeffErrNuMu', 'VeffErrNuTau']
	i=names.index(name)
	Column=DataArray[:,i]
	return Column

def RatiosWithErrors_Derp(Numerator,Denominator,nError,dError):
	Ratio=Numerator/Denominator
	RatioErrHi=(Numerator+nError)/(Denominator-dError) - Ratio
	RatioErrLo=Ratio - (Numerator-nError)/(Denominator+dError)
	return Ratio, RatioErrHi, RatioErrLo

def RatiosWithErrors(Numerator,Denominator,nError,dError):
	Ratio=Numerator/Denominator
	RatioErr=Ratio*np.sqrt(np.square(nError/Numerator) + np.square(dError/Denominator))
	return Ratio, RatioErr

def GetAllColumns(DataArray):
	E=GetColumn(DataArray,'E')
	VeffAve=GetColumn(DataArray,'VeffAve')
	VeffNuE=GetColumn(DataArray,'VeffNuE')
	VeffNuMu=GetColumn(DataArray,'VeffNuMu')
	VeffNuTau=GetColumn(DataArray,'VeffNuTau')
	VeffErrAve=GetColumn(DataArray,'VeffErrAve')
	VeffErrNuE=GetColumn(DataArray,'VeffErrNuE')
	VeffErrNuMu=GetColumn(DataArray,'VeffErrNuMu')
	VeffErrNuTau=GetColumn(DataArray,'VeffErrNuTau')
	return E, VeffAve, VeffNuE, VeffNuMu, VeffNuTau, VeffErrAve, VeffErrNuE, VeffErrNuMu, VeffErrNuTau

def SaveCSV(DataArray,outfn):
	names = 'E, VeffAve, VeffNuE, VeffNuMu, VeffNuTau, VeffErrAve, VeffErrNuE, VeffErrNuMu, VeffErrNuTau'
	np.savetxt(outfn,DataArray,delimiter=",",header=names)
