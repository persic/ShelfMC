import os,sys
import ROOT
import array

def getKey(nm):
    if (nm=="REFLECT_RATE"):
        return "R"
    elif (nm=="ATTEN_UP"):
        return "attUp"
    elif (nm=="ICETHICK"):
        return "iceThick"
    elif (nm=="FIRNfactor"):
        return "firnFactor"
    elif (nm=="reflect_events"):
        return "nReflect"
    elif (nm=="reflect_eventsweighted"):
        return "wtReflect"
    elif (nm=="direct_events"):
        return "nDirect"
    elif (nm=="direct_eventsweighted"):
        return "wtDirect"
    elif (nm=="combo_events"):
        return "nCombo"
    elif (nm=="combo_eventsweighted"):
        return "wtCombo"
    elif (nm=="nue"):
        return "nNuEtrg"
    elif (nm=="nuecounts"):
        return "nNuEthrown"
    elif (nm=="nueweighted"):
        return "wtNuE"
    elif (nm=="Veff_nue"):
        return "VeffNuE"
    elif (nm=="numu"):
        return "nNuMutrg"
    elif (nm=="numucounts"):
        return "nNuMuthrown"
    elif (nm=="numuweighted"):
        return "wtNuMu"
    elif (nm=="Veff_numu"):
        return "VeffNuMu"
    elif (nm=="nutau"):
        return "nNuTautrg"
    elif (nm=="nutaucounts"):
        return "nNuTauthrown"
    elif (nm=="nutauweighted"):
        return "wtNuTau"
    elif (nm=="Veff_nutau"):
        return "VeffNuTau"
    elif (nm=="Veff"):
        return "VeffAve"
    elif (nm=="integ_weight"):
        return "wtTot"
    elif (nm=="NNU"):
        return "nTotThrown"
    elif (nm=="AT"):
        return "nTrgAnts"
    elif (nm=="ST"):
        return "nStnAnts"
    elif (nm=="NSIGMA"):
        return "nSigmaThresh"
    else:
        return nm

def parseNVline(l, params):
    # extra space after NV=
    x = l.split("NV=")[-1].strip()
    x = x.split(' ')
    params["vthresh"] = float(x[0])
    x = x[1:]
    for y in x:
        kv = y.split('=')
        params[getKey(kv[0])] = float(kv[1])

def parseRefDirLine(l, params):
    x = l.split("=")
    y = [ a.split(" ") for a in x ]
    z = [ a for s in y for a in s ] # flatten
    z = filter(len, z) # get rid of empty strings
    q = [ a.split(",") for a in z ]
    for i in range(0, len(q), 2):
        params[getKey(q[i][0])] = int(q[i+1][0])   # N events
        params[getKey(q[i][0]+q[i][1])] = float(q[i+1][1]) # weight sum

def parseFlavorLine(l, params):
    x = l.split("=")
    y = [ a.split(" ") for a in x ]
    z = [ a for s in y for a in s ] # flatten
    z = filter(len, z) # get rid of empty strings
    q = [ a.split(",") for a in z ]
    z = [ a for s in q for a in s ] # flatten (again)
    z = filter(len, z) # get rid of empty strings
    for i in range(0, len(z), 6):
        nk,ck = z[i+0].split("/")
        nv,cv = z[i+3].split("/")
        ck = nk+ck # specify which nu count
        params[getKey(nk)]     = int(nv)
        params[getKey(ck)]     = int(cv)
        k = nk+z[i+1]
        params[getKey(k)]      = float(z[i+4])
        params[getKey(z[i+2])] = float(z[i+5])

def parseAveVeffLine(l, params):
    l = l.replace("= ","=")
    x = l.split(" ")
    x = filter(len, x) # get rid of empty strings
    x = filter(lambda y: '=' in y, x) # remove ones that don't have =
    for a in x:
        k,v = a.split("=")
        if (k == "AT/ST"):
            tk, nk = k.split("/")
            tv, nv = v.split("/")
            params[getKey(tk)] = int(tv)
            params[getKey(nk)] = int(nv)
        elif (k=="array"):
            c,r = v.split("x")
            params["ncols"] = int(c)
            params["nrows"] = int(r)
        elif ( (k=="NNU") or (k=="seed") ):
            params[getKey(k)] = int(v)
        else:
            params[getKey(k)] = float(v)

def parseEnergyLine(l, params):
    params["energy"] = float(l.split(" ")[0])


def parseVeffTextOutFile(infn):
    
    params = {}
    
    with open(infn) as inf:
        ignore = True
        for l in inf:
            l = l.strip()
            if (l.startswith("NV=")):
                # first line we want
                ignore = False

            if (ignore==False):
                if (l.startswith("NV=")):
                    parseNVline(l, params)
                elif (l.startswith("reflect")):
                    parseRefDirLine(l, params)
                elif (l.startswith("nue")):
                    parseFlavorLine(l, params)
                elif (l.startswith("Veff")):
                    parseAveVeffLine(l, params)
                elif (l.endswith("km3sr")):
                    parseEnergyLine(l, params)
                else:
                    print "Ingoring line [{0}]".format(l)

            if (l.endswith("km3sr")):
                # don't need any more lines
                ignore = True

    return params

def parseVeffTextOut(indir, outfn):
    
    outf = ROOT.TFile(outfn, "recreate")
    outf.cd()
    nt = ROOT.TTree("nt","shelfmc Veff parameters")
    
    brs = {}
    
    makeBranches = True
    for path,dirs,files in os.walk(indir):
        for pinfn in files:
            infn = path+"/"+pinfn
            print "Parsing [{0}]".format(infn)
            params = parseVeffTextOutFile(infn)
            #print params
            if (len(params)>0):
                if (makeBranches):
                    for k,v in params.iteritems():
                        tp = 'f'
                        tb = '/F'
                        if (type(v)==int):
                            tp = 'i'
                            tb = '/I'
                        x = array.array(tp, [0])
                        brs[k] = x
                        nt.Branch(k, x, k+tb)
                    makeBranches = False
                
                for k,v in params.iteritems():
                    brs[k][0] = v
                nt.Fill()

    outf.Write()
    outf.Close()
    
    print "Wrote file [{0}]".format(outfn)

    
if __name__=="__main__":
    
    if (len(sys.argv)<3):
        print "Usage: python parseVeffTextOut [input directory] [output file]"
        sys.exit(1)

    parseVeffTextOut(sys.argv[1], sys.argv[2])

