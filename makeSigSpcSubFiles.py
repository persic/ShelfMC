import os,sys
import numpy
import subprocess
import time

quetag         = "AGKZS" # must start with letter
basedir        = "/pub/cpersich/ShelfMCGit/ShelfMC/ARIANNA_GKZ_Shadow"
runVeffShelfMC = "/pub/cpersich/ShelfMCGit/ShelfMC/runVeffShelfMC.py"

MaxEvntPerJob=int(1e5)
NNUTotal=int(1e7)
startSeed=41

inrefn = "input_reference.txt"

def checkDirExists(d):
    if not os.path.exists(d):
        try:
            os.makedirs(d)
        except OSError as ex:
            if ex.errno != errno.EEXIST:
                raise
SpecName = {
0:"E-2",
1:"GKZ"
}

SPECTRUM= 1 #Use Spectrum if 1, single energy if 0
WIDESPECTRUM = 1 #Exp 16-21.5 if 1, 17-20 if 0
GKZ = 1 #Use GKZ spectrum if 1, E-2 if 0
ICETHICK=575 #575 for Moore's Bay, 2700 for SP
ATTEN_UP=500 #Average Atten Length (or something like that)
refl  = [ 0.9 ] #reflectivity
shadowing = 1
ATGap = 1000 #distance between stations in meters
ST4R = 3.0 #Radius in meters between center station and antenna
FREQ_LOW = 50 #low frequency of LPDA Response MHz
FREQ_HIGH = 1000 #high frequency of LPDA Response MHz
SIGNAL_FLUCT = 1 #1=add noise fluctuation to signal or 0=do not

# nrows, ncols, n_ant_perST, n_ant_trigger, nsigma
stations = [ (1,1,4,2,4),  # single Stn
             #(3,5,4,2,4),  # HRA-7
             #(7,7,6,3,3),  # 49 Stn
             #(36,36,4,2,4), # full ARIANNA
             #(31,31,8,3,5) # 960 stn ARIANNA, 5sigma
             ]

ilines = []
with open(inrefn) as inref:
    for l in inref:
        ilines.append(l)

j = 1

for r in refl:
    for stnpars in stations:
        nrows = stnpars[0]
        ncols = stnpars[1]
        nants = stnpars[2]
        natrg = stnpars[3]
        nsigs = stnpars[4]

        if NNUTotal <= MaxEvntPerJob:
            nsubjobs=1
        elif NNUTotal > MaxEvntPerJob:
            if NNUTotal%MaxEvntPerJob == 0:
                    nsubjobs=NNUTotal/MaxEvntPerJob
            else:
                    nsubjobs=NNUTotal/MaxEvntPerJob + 1

        for i in range(1,nsubjobs+1):
                seed=startSeed+j
                if i*MaxEvntPerJob <= NNUTotal:
                    n=MaxEvntPerJob
                elif i*MaxEvntPerJob > NNUTotal:
                    n=NNUTotal-(i-1)*MaxEvntPerJob
                outdir = "SigSpace_Spec{spec}_r{r:0.1f}_stn{row}-{col}_sig{sig}_"\
                    "trg{mj}-{ant}_subjob{i}of{nsubjobs}".format(spec=SpecName[GKZ], r=r,
                                           row=nrows, col=ncols,
                                           sig=nsigs, mj=natrg, ant=nants, i=i, nsubjobs=nsubjobs)

                foutdir = os.path.join(basedir, "out", outdir)

                checkDirExists(foutdir)

                infn = os.path.join(foutdir, "input.txt")
                with open(infn,"w") as inf:
                    for l in ilines:
                        if   ("#NNU" in l):
                            l = "{0:d} #NNU, setting to 1 for unique neutrino\n"\
                                .format(n)
                        if   ("#seed" in l):
                            l = "{0:d}      #seed Seed for Rand3\n"\
                                .format(seed)
                        elif ("#EXPONENT" in l):
                            l = "{0:0.1f}    #EXPONENT, !should be exclusive with "\
                                "SPECTRUM\n".format(0)
                        elif ("#ATGap" in l):
		            l = "{0:0.1f}    #ATGap, distance between stations in meters\n".format(ATGap)
                        elif ("#N_Ant_perST" in l):
                            l = "{0:d} #N_Ant_perST, not to be confused with "\
                                "ST_TYPE above\n".format(nants)
                        elif ("#N_Ant_Trigger" in l):
                            l = "{0:d} #N_Ant_Trigger, this is the minimum number "\
                                "of AT to trigger\n".format(natrg)
                        elif ("#NROWS" in l):
                            l = "{0:d} #NROWS 12 initially, set to 3 for "\
                                "HEXAGONAL\n".format(nrows)
                        elif ("#NCOLS" in l):
                            l = "{0:d} #NCOLS 12 initially, set to 5 for "\
                                "HEXAGONAL\n".format(ncols)
                        elif ("#NSIGMA" in l):
                            l = "{0:d} #NSIGMA, threshold of trigger\n"\
                                .format(nsigs)
                        elif ("#SPECTRUM" in l):
                            l = "{0:d} #SPECTRUM, use spectrum, ! was 1 initially!\n"\
                                .format(SPECTRUM)
                        elif ("#ICETHICK" in l):
                            l = "{0:0.1f}  #ICETHICK, thickness of ice including firn, "\
                                "575m at Moore's Bay\n".format(ICETHICK)
                        elif ("#ATTEN_UP" in l):
                            l = "{0:0.1f}, this is the conjuction of the plot attenlength_up "\
                                "and attlength_down when setting REFLECT_RATE=0.5(3dB)\n"\
                                .format(ATTEN_UP)
                        elif ("#REFLECT_RATE" in l):
                            l = "{0:0.1f}    #REFLECT_RATE,power reflection rate "\
                                "at the ice bottom\n".format(r)
                        elif   ("#SHADOWING" in l):
                            l = "{0:d}      #SHADOWING\n"\
                                .format(shadowing)
                        elif ("#GKZ" in l):
                            l = "{0:d} #GZK, 1 means using GZK flux, 0 means E-2 flux\n"\
                                .format(GKZ)
                        elif ("#WIDESPECTRUM" in l):
                            l = "{0:d} #WIDESPECTRUM, use 10^16 eV to 10^21.5 eV as the energy "\
                                "spectrum, otherwise use 17-20\n".format(WIDESPECTRUM)
                        elif ("#HEXAGONAL" in l):
                            hexag = 0
                            if ( (nrows==3) and (ncols==5) ):
                                hexag = 1
                            l = "{0:d} #HEXAGONAL\n".format(hexag)
                        elif ("#ST4_R" in l):
                            l = "{0:0.1f}     #ST4_R radius in meters between center "\
                                "of station and antenna\n".format(ST4R)
                        elif   ("#SIGNAL_FLUCT" in l):
                            l = "{0:d}      #SIGNAL_FLUCT 1=add noise fluctuation to signal or 0=do not\n"\
                                .format(SIGNAL_FLUCT)
                        elif ("#FREQ_LOW" in l):
                            l = "{0:0.1f}     #FREQ_LOW low frequency of LPDA Response MHz\n"\
                                .format(FREQ_LOW)
                        elif ("#FREQ_HIGH" in l):
                            l = "{0:0.1f}     #FREQ_HIGH high frequency of LPDA Response MHz\n"\
                                .format(FREQ_HIGH)
                        # else just copy the line without modification
                        inf.write(l)

                jobdir = os.path.join(basedir, "jobs")
                jobfn = os.path.join(jobdir, "shelfmcJob.{0}".format(j))
                checkDirExists(jobdir)
                with open(jobfn,"w") as jf:
                    jf.write("{foutdir} {outprfx}".format(
                            foutdir=foutdir, outprfx=outdir))
                j+=1


# make the queue file
logd = os.path.join(basedir, "logs")
checkDirExists(logd)
with open( os.path.join(basedir,
                        "queueSigSpc{quetag}.sh".format(quetag=quetag)),
           "w" ) as qf:
    qf.write(
        "#!/bin/bash\n"
        "# use an array job\n"
        "#$ -t 1-{njobs}\n"
        "#$ -N {quetag}Veff\n"
        "#$ -q grb,grb64,pub64,free64\n"
        "#$ -V\n"
        "#$ -o {logdir}\n"
        "#$ -j y\n"
        "#$ -ckpt restart\n"
        "\n"
        "ulimit -c 0\n"
        "\n"
        "python {runVeffShelfMC} "\
            "{jobdir}/shelfmcJob.${{SGE_TASK_ID}}\n"
        "\n"
        .format(njobs=j-1, quetag=quetag, logdir=logd, 
                runVeffShelfMC=runVeffShelfMC, jobdir=jobdir))


