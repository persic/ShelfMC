import sys, os
import xml.etree.ElementTree as ET
import numpy as np

def indent(elem, level=0):
    indentStr = "    "
    i = "\n" + level*indentStr
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + indentStr
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

#Apparently this is the wrong convention
#def CartFromPolar(r, theta, phi):
#    x=r*np.sin(theta)*np.cos(phi)
#    y=r*np.sin(theta)*np.sin(phi)
#    z=r*np.cos(theta)
#    return x,y,z
#
#def PolarFromCart(x,y,z):
#    r = np.linalg.norm([x,y,z])
#    theta = np.arccos(z/r)
#    phi = np.arctan2(y,x)
#    return r,theta,phi

#Antenna positions from ARASim uses the convention -pi/2<Theta<pi/2
def CartFromPolarPosNu(r, theta, phi):
    x=r*np.cos(theta)*np.cos(phi)
    y=r*np.cos(theta)*np.sin(phi)
    z=r*np.sin(theta)
    return x,y,z

def PolarFromCartPosNu(x,y,z):
    r = np.linalg.norm([x,y,z])
    theta = np.arcsin(z/r)
    phi = np.arctan2(y,x)
    return r,theta,phi

def RandomNormal(D): #D is the number of dimensions
    Vect = np.random.normal(size=D)
    return Vect/np.linalg.norm(Vect)

def RandomVectorSphere(D): #D is the number of dimensions
    UVector = RandomNormal(D)
    R = np.random.rand()
    return UVector * (R**(1.0/D))

def RandomPhiTheta():
    phi = 2*np.pi*np.random.rand()
    theta = np.arccos(2*np.random.rand()-1)
    return phi, theta

class EventProps:
    FlavorMap = {
    'nue':1,
    'numu':2,
    'nutau':3
    }
    CurrentMap = {
    'cc':0,
    'nc':1
    }
    AntiMap = {
    'nu':0,
    'nubar':1
    }

    def __init__(self,Num=None,Exponent=None,Flavor=None,Inelasticity=None,\
    Current=None,X=None,Y=None,Z=None,RPos=None,ThetaPos=None,PhiPos=None,
    ThetaDir=None,PhiDir=None,NuNuBar=None):
        self._Num=Num
        self._Exponent = Exponent
        self._Flavor = Flavor
        self._Inelasticity = Inelasticity
        self._Current = Current
        self._X = X
        self._Y = Y
        self._Z = Z
        self._RPos = RPos
        self._ThetaPos = ThetaPos
        self._PhiPos = PhiPos
        self._ThetaDir = ThetaDir
        self._PhiDir = PhiDir
        self._NuNuBar = NuNuBar
        self._XYZUpToDate = (self._X!=None and self._Y!=None and self._Z!=None)
        self._RThetaPhiUpToDate = (self._RPos!=None and self._ThetaPos!=None and self._PhiPos!=None)

    def UpdateXYZ(self):
        if (not self._XYZUpToDate) and self._RThetaPhiUpToDate:
            self._X, self._Y, self._Z = CartFromPolarPosNu(self._RPos,self._ThetaPos,self._PhiPos)
            self._XYZUpToDate = True
        elif (not self._XYZUpToDate) and (not self._RThetaPhiUpToDate):
            print('WARNING: UpdateXYZ Coordinates Never Set')

    def UpdateRThetaPhi(self):
        if (not self._RThetaPhiUpToDate) and self._XYZUpToDate:
            self._RPos, self._ThetaPos, self._PhiPos, = PolarFromCartPosNu(self._RPos,self._ThetaPos,self._PhiPos)
            self._RThetaPhiUpToDate = True
        elif (not self._XYZUpToDate) and (not self._RThetaPhiUpToDate):
            print('WARNING: UpdateRThetaPhi Coordinates Never Set')

    def ShiftOrigin(self,Origin):
        if self._XYZUpToDate == False:
            self.UpdateXYZ()
        self._X -= Origin[0]
        self._Y -= Origin[1]
        self._Z -= Origin[2]
        self._XYZUpToDate = True
        self._RThetaPhiUpToDate = False

    def GetNum(self):
        return self._Num
    def GetExponent(self):
        return self._Exponent
    def GetFlavorString(self):
        return self._Flavor
    def GetFlavorInt(self):
        if self._Flavor in self.FlavorMap:
            return self.FlavorMap[self._Flavor]
        else:
            print('Invalid Flavor String in Event {}'.format(self._Num))
            return None
    def GetInelasticity(self):
        return self._Inelasticity
    def GetCurrentString(self):
        return self._Current
    def GetCurrentInt(self):
        if self._Current in self.CurrentMap:
            return self.CurrentMap[self._Current]
        else:
            print('Invalid Current String in Event {}'.format(self._Num))
            return None
    def GetXYZ(self):
        if self._XYZUpToDate:
            return self._X, self._Y, self._Z
        elif self._RThetaPhiUpToDate:
            self._X, self._Y, self._Z = CartFromPolarPosNu(self._RPos,self._ThetaPos,self._PhiPos)
            self._XYZUpToDate=True
            return self._X, self._Y, self._Z
        else:
            print('Coordinates Not Set')
            return None
    def GetRThetaPhi(self):
        if self._RThetaPhiUpToDate:
            return self._RPos, self._ThetaPos, self._PhiPos
        elif self._XYZUpToDate:
            self._RPos, self._ThetaPos, self._PhiPos = PolarFromCartPosNu(self._X,self._Y,self._Z)
            self._RThetaPhiUpToDate = True
            return self._RPos, self._ThetaPos, self._PhiPos
        else:
            print('Coordinates Not Set')
            return None
    def GetThetaDir(self):
        return self._ThetaDir
    def GetPhiDir(self):
        return self._PhiDir
    def GetNuNuBarSting(self):
        return self._NuNuBar
    def GetNuNuBarInt(self):
        if self._NuNuBar in self.AntiMap:
            return self.AntiMap[self._NuNuBar]
        else:
            print('Invalid nu/nubar String in Event {}'.format(self._Num))
            return None

    def SetNum(self,Num):
        self._Num = int(Num)
    def SetExponent(self,Exponent):
        self._Exponent = Exponent
    def SetFlavorString(self,Flavor):
        if Flavor in self.FlavorMap:
            self._Flavor = Flavor
        else:
            print('Invalid flavor string')
    def SetFlavorInt(self,FlavorInt):
        InMap = False
        for FlavorString, Int in self.FlavorMap.iteritems():
            if Int == FlavorInt:
                InMap = True
                self._Flavor = FlavorString
        if InMap==False:
            print('Invalid flavor id')
    def SetCurrentString(self,Current):
        if Current in self.CurrentMap:
            self._Current = Current
        else:
            print('Invalid current string')
    def SetCurrentInt(self,CurrentInt):
        InMap = False
        for CurrentString, Int in self.CurrentMap.iteritems():
            if Int == CurrentInt:
                InMap = True
                self._Current = CurrentString
        if InMap==False:
            print('Invalid current id')
    def SetInelasticity(self,Inelasticity):
        self._Inelasticity = Inelasticity
    def SetX(self,X):
        self._X = X
        self._XYZUpToDate = (self._X!=None and self._Y!=None and self._Z!=None)

    def SetY(self,Y):
        self._Y = Y
        self._XYZUpToDate = (self._X!=None and self._Y!=None and self._Z!=None)

    def SetZ(self,Z):
        self._Z = Z
        self._XYZUpToDate = (self._X!=None and self._Y!=None and self._Z!=None)

    def SetRPos(self,RPos):
        self._RPos = RPos
        self._RThetaPhiUpToDate = (self._RPos!=None and self._ThetaPos!=None and self._PhiPos!=None)
    def SetThetaPos(self,ThetaPos):
        self._ThetaPos = ThetaPos
        self._RThetaPhiUpToDate = (self._RPos!=None and self._ThetaPos!=None and self._PhiPos!=None)
    def SetPhiPos(self,PhiPos):
        self._PhiPos = PhiPos
        self._RThetaPhiUpToDate = (self._RPos!=None and self._ThetaPos!=None and self._PhiPos!=None)
    def SetThetaDir(self,ThetaDir):
        self._ThetaDir = ThetaDir
    def SetPhiDir(self,PhiDir):
        self._PhiDir = PhiDir
    def SetNuNuBar(self,NuNuBar):
        self._NuNuBar = NuNuBar
    def SetNuNuBarString(self,NuNuBar):
        if NuNuBar in self.AntiMap:
            self._NuNuBar = NuNuBar
        else:
            print('Invalid nu/nubar string')
    def SetNuNuBarInt(self,NuNuBarInt):
        InMap = False
        for NuNuBarString, Int in self.AntiMap.iteritems():
            if Int == NuNuBarInt:
                InMap = True
                self._NuNuBar = NuNuBarString
        if InMap==False:
            print('Invalid current id')

def ReadEventListXMLToListObject(infn):
    tree = ET.parse(infn)
    root = tree.getroot()
    Events=[]
    for List in root.iter('List'):
        for i, Event in enumerate(List.iter('Event')):

            E = EventProps()
            EvID = int(Event.find('EvID').text)
            E.SetNum(EvID)

            Exponent = float(Event.find('Exponent').text)
            E.SetExponent(Exponent)

            Flavor = Event.find('Flavor').text
            E.SetFlavorString(Flavor)

            Current = Event.find('Current').text
            E.SetCurrentString(Current)

            Inelasticity = float(Event.find('Inelasticity').text)
            E.SetInelasticity(Inelasticity)

            X = float(Event.find('X').text)
            E.SetX(X)

            Y = float(Event.find('Y').text)
            E.SetY(Y)

            Z = float(Event.find('Z').text)
            E.SetZ(Z)

            ThetaDir = float(Event.find('ThetaDir').text)
            E.SetThetaDir(ThetaDir)

            PhiDir = float(Event.find('PhiDir').text)
            E.SetPhiDir(PhiDir)

            Events.append(E)

    return Events

def ReadEventListARASimToListObject(infn):
    Events=[]
    inArray = np.genfromtxt(infn,comments='//',skip_header=5)

    if len(inArray.shape)==1:
        inArray = [inArray,]

    for (Num, FlavorInt, NuNuBar, Exponent, CurrentInt, RPos, ThetaPos, PhiPos, ThetaDir, PhiDir, Inelasticity) in inArray:

        E = EventProps()

        E.SetNum(Num)
        E.SetExponent(Exponent)
        E.SetFlavorInt(FlavorInt)
        E.SetCurrentInt(CurrentInt)
        E.SetNuNuBarInt(NuNuBar)
        E.SetRPos(RPos)
        E.SetThetaPos(ThetaPos)
        E.SetPhiPos(PhiPos)
        E.SetThetaDir(ThetaDir)
        E.SetPhiDir(PhiDir)
        E.SetInelasticity(Inelasticity)

        Events.append(E)
    return Events

def WriteARASimList(outfn,Events,Origin=[0,0,2700.]):
    NumEvents = len(Events)

    header = '//version\nVERSION=0.1\nEVENT_NUM={}\n//evid nuflavorint nu_nubar pnu currentint posnu_r posnu_theta posnu_phi nnu_theta nnu_phi elast_y\n//use -999 as a flag of random input'.format(NumEvents)

    typeFormat = ['%i','%i','%i','%f','%i','%f','%f','%f','%f','%f','%f']

    OutArray = []
    for E in Events:
        E.ShiftOrigin(Origin)
        Num = E.GetNum()
        FlavorInt = E.GetFlavorInt()
        if E.GetNuNuBarSting() == None:
            NuNuBar = 0
        else:
            NuNuBar = E.GetNuNuBarInt()
        Exponent = E.GetExponent()
        CurrentInt = E.GetCurrentInt()
        RPos, ThetaPos, PhiPos = E.GetRThetaPhi()
        ThetaDir = E.GetThetaDir()
        PhiDir = E.GetPhiDir()
        Elasticity = E.GetInelasticity()

        OutArray.append((Num,FlavorInt,NuNuBar,Exponent,CurrentInt,RPos,ThetaPos,PhiPos,ThetaDir,PhiDir,Elasticity))
    OutArray = np.array(OutArray)
    np.savetxt(outfn,OutArray,header=header,fmt=typeFormat,delimiter=' ')

def ARASimListArray(Events,Origin=[0,0,0.0]):
    NumEvents = len(Events)

    header = '//version\nVERSION=0.1\nEVENT_NUM={}\n//evid nuflavorint nu_nubar pnu currentint posnu_r posnu_theta posnu_phi nnu_theta nnu_phi elast_y\n//use -999 as a flag of random input'.format(NumEvents)

    typeFormat = ['%i','%i','%i','%f','%i','%f','%f','%f','%f','%f','%f']

    OutArray = []
    for E in Events:
        E.ShiftOrigin(Origin)
        Num = E.GetNum()
        FlavorInt = E.GetFlavorInt()
        if E.GetNuNuBarString() == None:
            NuNuBar = 0
        else:
            NuNuBar = E.GetNuNuBarInt()
        Exponent = E.GetExponent()
        CurrentInt = E.GetCurrentInt()
        RPos, ThetaPos, PhiPos = E.GetRThetaPhi()
        ThetaDir = E.GetThetaDir()
        PhiDir = E.GetPhiDir()
        Elasticity = E.GetInelasticity()

        OutArray.append((Num,FlavorInt,NuNuBar,Exponent,CurrentInt,RPos,ThetaPos,PhiPos,ThetaDir,PhiDir,Elasticity))
    OutArray = np.array(OutArray)
    return OutArray

def WriteShelfMCXMLList(outFN,Events,Origin=[0,0,-2700]):
    EventList = ET.Element('EventList')
    List = ET.SubElement(EventList,'List')

    for E in Events:
        E.ShiftOrigin(Origin)

        EventElemnt = ET.SubElement(List,'Event')

        EvID = E.GetNum()
        EvIDElement = ET.SubElement(EventElemnt,'EvID')
        EvIDElement.text = '{:d}'.format(int(EvID))

        Exponent = E.GetExponent()
        ExponentElement = ET.SubElement(EventElemnt,'Exponent')
        ExponentElement.text = '{}'.format(Exponent)

        Flavor = E.GetFlavorString()
        FlavorElement = ET.SubElement(EventElemnt,'Flavor')
        FlavorElement.text = Flavor

        Inelasticity = E.GetInelasticity()
        InelasticityElement = ET.SubElement(EventElemnt,'Inelasticity')
        InelasticityElement.text = '{}'.format(Inelasticity)

        Current = E.GetCurrentString()
        CurrentElement = ET.SubElement(EventElemnt,'Current')
        CurrentElement.text = Current

        X, Y, Z = E.GetXYZ()
        XElement = ET.SubElement(EventElemnt,'X')
        XElement.text = '{}'.format(X)

        YElement = ET.SubElement(EventElemnt,'Y')
        YElement.text = '{}'.format(Y)

        ZElement = ET.SubElement(EventElemnt,'Z')
        ZElement.text = '{}'.format(Z)

        ThetaDir = E.GetThetaDir()
        ThetaDirElement = ET.SubElement(EventElemnt,'ThetaDir')
        ThetaDirElement.text = '{}'.format(ThetaDir)

        PhiDir = E.GetPhiDir()
        PhiDirElement = ET.SubElement(EventElemnt,'PhiDir')
        PhiDirElement.text = '{}'.format(PhiDir)

    indent(EventList)
    Tree = ET.ElementTree(element=EventList)
    Tree.write(outFN)

def WriteNuRadioMCimList(outfn,Events,Origin=[0,0,0]):
    NumEvents = len(Events)
    Version = "0.1"
    FlavorMap = {
    'nue':'e',
    'numu':'mu',
    'nutau':'tau'
    }
    header = 'VERSION={}'.format(Version)

    typeFormat = ['%i','%s','%s','%e','%s','%f','%f','%f','%f','%f','%f']

    OutArray = []
    for E in Events:
        E.ShiftOrigin(Origin)
        Num = E.GetNum()
        FlavorString = FlavorMap[E.GetFlavorString()]
        if E.GetNuNuBarSting() == None:
            NuNuBar = 'nu'
        else:
            NuNuBar = E.GetNuNuBarSting()
        Exponent = E.GetExponent()
        E = 10**Exponent
        CurrentString = E.GetCurrentString()
        X, Y, Z = E.GetXYZ()
        ThetaDir = E.GetThetaDir()
        PhiDir = E.GetPhiDir()
        Elasticity = E.GetInelasticity()

        OutArray.append((Num,FlavorInt,NuNuBar,Exponent,CurrentInt,RPos,ThetaPos,PhiPos,ThetaDir,PhiDir,Elasticity))
    OutArray = np.array(OutArray)
    np.savetxt(outfn,OutArray,header=header,fmt=typeFormat)

def makeRandomEventList(N, Exp, RhoMax = 4000.0, IceThick=2700.0, nunubar = 1):
    """Creates Random neutrino vertices in ShelfMC's Coordinate system, with the origin at the bottom of the ice. """
    Events = []

    for i in xrange(N):
        x,y = RhoMax*RandomVectorSphere(2)
        z = IceThick*np.random.rand()
        #print('[x,y,z] = [{},{},{}]'.format(x,y,z))
        gamma = np.random.rand()
        flavor=np.random.randint(1,4)
        current=np.random.randint(2)
        phi, theta = RandomPhiTheta()

        E = EventProps(
        Num=i,
        Exponent=Exp,
        Inelasticity = gamma,
        NuNuBar=nunubar,
        X=x,Y=y,Z=z,
        PhiDir=phi,ThetaDir=theta
        )

        E.SetFlavorInt(flavor)
        E.SetCurrentInt(current)

        Events.append(E)

    return Events

if __name__ == '__main__':
    infn = sys.argv[1]
    out1 = sys.argv[2]
    out2 = sys.argv[3]

    print('reading in {}'.format(infn))
    Events = ReadEventListXMLToListObject(infn)

    print('converting to ARA format and writing to {}'.format(out1))
    WriteARASimList(out1,Events)

    print('reading in {}'.format(out1))
    Events = ReadEventListARASimToListObject(out1)

    print('converting to ShelfMC format and writing to {}'.format(out2))
    WriteShelfMCXMLList(out2,Events)
