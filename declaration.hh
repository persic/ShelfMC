// CJR 2015-07-15: add gain file as input.txt paramter
# include "tinyxml2.h"
const char * StnGeoFN;

//constant parameters
const double C = 3.e8; //speed of light m/s
const double R_earth = 6357390.; //m
const double INCH2M = 0.0254; // 1 inch=0.0254m
const double SEC2YR = 1. / (365.*24.*3600.); // convert second to year
const double YR2SEC = 365 * 24.*3600; //year to seconds
const double PI = 3.14159265;
const double ALOG2 = 0.693147;
const double RAD2DEG = 57.2957795;
const double DEG2RAD = 0.0174532925;
const double NICE = 1.78;
const double X0ICE = 0.403; //parameter from Jaime's paper
const double ECICE = 63.7; //critical energy in ice(MeV)
const double AEX_ICE = 1.; //efficiency for producing charge asymmetry relative to ice. 1 by definition
const double DensityAIR = 1.25; //kg/m^3
const double DensityICE = 917; //kg/m^3
const double DensityWATER = 1000.; //kg/m^3
const double DensityCRUST = 2900.; //kg/m^3 //SB 11/18/11 change to 3500 to better reflect mantle
const double AMU = 1.66E-27; //kg
double JAIME_FACTOR = 1.; //factor to multiply Jaime's parameterization for error analysis
const double SIGMA_FACTOR = 1.;
const double LIVETIME = 180.*24.*3600.; //half year in seconds//added factor of 2 on 9/13

int seed;
double FIRNDEPTH;//m
double NSIGMA;//threshold for trigger
double TNOISE;//Noise Temperature
double EXPONENT;//energy exponent
int CONST_ATTENLENGTH;//0 means not constant attenlength
int ATTEN_EXP = 1; //use the data from Minna Bluff measurement
int ATTEN_FREQ = 0; //let attenuation length be a function of frequency when f<122MHz
double ATTEN_UP;
double ATTEN_DOWN;
int SIGNAL_FLUCT = 0; //1=add noise fluctuation to signal or 0=do not
double ICETHICK; //575 for Moore's Bay
int FIRN;//0 means only uniform ice, 1 means ice and firn 2 uniform layers.
int SCATTER;
double SCATTER_WIDTH;
int SPECTRUM;
int WIDESPECTRUM;
int DIPOLE;
int N_Ant_perST = 0; //KD now input as part of input file
int N_Ant_Trigger; //KD new majority logic variable
int n_chan_perST;//the total number of channels each station
double  NFIRN;//the initial refraction index at the surface of the ice
int NNU;//number of neutrinos isotropically reach the ice block
double ATTEN_FACTOR;
double REFLECT_RATE;
int GZK;//whether use GZK flux of E-2 flux
int FANFLUX;
double  Z;//downward angle of 4 side antennas when ST_TYPE=1
int DEPTH_DEPENDENT_N; //shifted to steering file
//int DEPTH_DEPENDENT_N=0;//0 means uniform firn, 1 means n of firn is a function of depth
int LPM = 1;
int SECKEL = 0;
double n1 = 1.78; //ice refraction index at the interaction position

int TAUREGENERATION = 1; //KD: if 1=tau regeneration effect, if 0=original
int PLANEWAVE = 1; //KD: if 1: plane wave calculations is on //check code for multistation
//int SHADOWING = 0; //if FIRN is 0, it automatically doesn't apply
int SHADOWING; //now shifted to steering file 11/17/11
int HEXAGONAL = 0; //note changes in row,col(input file) and ATGAPs(function file) (500m, 866m), works for 7 stations
double HRAfactor = 1; //sets the spacing scale wrt to 1km
double FIRNfactor = 1;

//parameters of antennas
//const double HeightAT=37.62*INCH2M;//height of antenna in m
//const double WidthAT=37.62*INCH2M;//width of antenna in m
//const double LengthAT=22.63*INCH2M;//length of antenna in m
int NROWS;
int NCOLS;
double ATGap;//m
double ST_DISTANCE = 4.; //m
//double ST4_R=1./sin(22.5*DEG2RAD);
double ST4_R = 3.5;
//const double ATGap=200.;//m
double FREQ_LOW; //MHz
double FREQ_HIGH; //MHz;
double BW; //MHz Now Declaired in shelfmc_stripped.cc
double FREQ_BIN;
double EDGE;//m
//const double ICETHICK = 575.; //m //FW:624 //double ICETHICK=575.; //Moved to input.txt
const double ATDepth = 0.; //put the antennas at the mirror surface of the real ice surface taking the ice bottom as a mirror.


//parameters of seavey
const double heff_max = 0.62639; //m,effective height of seavey antenna
double heff_max_LPA = heff_max; //suppose the longest dipole in the LPA is 2 meters, 2*0.64m as the maximum effective height
int NPOINTS; // number of GPS positions we're picking from.
int NPOINTS_GAIN; // number of points within bandwidth that gain is measured

const double bwslice_center[4] = {265.e6, 435.e6, 650.e6, 980.e6 }; // center frequencies
const double bwslice_width[4] = {130.e6, 160.e6, 250.e6, 370.e6}; // 3 dB bandwidths, without overlap

//parameters of dipole from Amy
const double gain_dipole = 2.15;
const double Z0 = 377.; //ohms, intrinsic impedance of free space(a pure resistance)
const double Zr = 50.; //ohms,radation resistance

//signals
const double CONEWIDTH = 20.*DEG2RAD;
const int NFREQ = 95; //number of frequency bins 95//380
const int NVIEWANGLE = 7; //number of viewing angles to look at the signal,on Seckel's request
//const double FREQ_LOW=10;//MHz
//const double FREQ_HIGH=5000;//MHz Assigned in input.txt 
//const double FREQ_LOW = 100; //MHz
//const double FREQ_HIGH = 1000; //MHz;
//const double BW = FREQ_HIGH - FREQ_LOW; //MHz Now Declaired in shelfmc_stripped.cc
//const double FREQ_BIN = BW / NFREQ;
const double KBOLTZ = 1.38E-23; //Boltzman constant J/K

double pnu;//energy of neutrino beginning for 1E17 ev, a little bit lower than the energy in ANITA
double changle;//cherenkov angle


//KD 08/10
double gainv; //KD for gain dependency
//double flare_v; //KD for gain dependency
//double flare_h; //KD for gain dependency
//double GetGainquick(double flare_v, double flare_h, double gainv); //KD for gain dependency


double flare[4][NFREQ];
double gain[2][NFREQ];
double gainv_measured[100000];
double gainh_measured[100000];
double frequency_forgain_measured[100000];

//06/26 at guesthouse
double vmmhz1m_em; // V/m/MHz at 1m due to EM component of shower
double vmmhz1m_had; // V/m/MHz at 1m due to HAD component of shower

double freq[NFREQ]; //freq for each bin, MHz

double forseckel[NVIEWANGLE][NFREQ];// Per Seckel's request, get strength of signal across frequencies for different viewing angles.
double MINVIEWANGLE = acos(1 / NICE) - 0.3;
double MAXVIEWANGLE = acos(1 / NICE) + 0.3;
double viewangles[NVIEWANGLE] = {acos(1 / NICE),
                                 acos(1 / NICE) - 2.*DEG2RAD,
                                 acos(1 / NICE) - 5.*DEG2RAD,
                                 acos(1 / NICE) - 10.*DEG2RAD,
                                 acos(1 / NICE) - 20.*DEG2RAD,
                                 acos(1 / NICE) - 30.*DEG2RAD,
                                 acos(1 / NICE) - 40.*DEG2RAD
                                };

struct AntennaPlacement{
  int Type;
  double position[3];
  double n_boresight[3];
  double n_eplane[3];
};

double GetNoise(double);
double GetDecayLength(double EXPONENT); //CP 07/15
void GetBothLocation(int, int, double*, double*);
void GetATLocation(int, int, double*);
void GetMirrorATLocation(int, int, double*);
double GetGainV(double freq);
double GetGainH(double freq);
double GaintoHeight(double gain, double freq);

double GetHeff(int, double, double*, double*, double*, double*);

double GetN(double* posnu);
double GetRange(double firndepth); //KD added 9/23/10
void GetInteractionPoint(double* posnu);
void GetPosnuBin(double* posnu, int& posnu_iRow, int& posnu_iCol);
void GetAttenlength(double* posnu, double& attenlength_up, double& attenlength_down);
void GetEntryPoint(double theta_nu, double* nnu, double* posnu, double* entry);

void GetBouncePoint(double theta_nposnu2MirrorAT, double* nposnu2MirrorAT, double* posnu, double* bounce); //KD 06/29/11

void GetChord(double theta_nu, double& chord_inCRUST, double& chord2_inICE);
double GetWeight(double, double, double, double, double);
void Sigma(double pnu, double& sigma, double& int_len_kgm2);

void GetVmMHz(double, double, double, double*);
void GetVmMHz_freq(double*, double, double, double*); //make vmmhz_max a function of frequency since the attenuation length is a function of frequency
double GetVmMHz1m(double, double, double, double, double, double);
void TaperVmMHz(double viewangle, double deltheta_em, double deltheta_had, double emfrac, double hadfrac,
                double& vmmhz,
                double& vmmhz_em, //KD for energy tracking
                double& vmmhz_had //KD for energy tracking
               );
void GetSpread(double pnu,
               double emfrac,
               double hadfrac,
               double freq,
               double n_depth,
               double X0DEPTH,

               double& deltheta_em_max,
               double& deltheta_had_max,

               double& eshower_em,//KD
               double& eshower_had, //KD
               double& showerlength_em //KD
              );
double VmMHz1m(double viewangle, double freq, double emfrac, double hadfrac, double em_shower_length, double had_shower_length); //Seckel's scaling
double VmMHz(double vmmhz1m, double r);
double VmMHz_attenuated(double d_posnu2AT, double vmmhz, double attenlength);
void GetFlare(double freq, double* flare_tmp);
int GetBeamWidths(double flare[4][NFREQ], double gain[2][NFREQ], double freq[NFREQ]);
void ReadGains();
void ConvertEHtoLR(double e_component, double h_component, double& lcp_component, double& rcp_component);
void GetPolarization(double* nnu, double* posnu2AT, double* n_pol, double* n_Bfield);

//KD adapted from ANITA MC 08/25
void GetFresnel(double theta1, double theta2, double* original_nsignal, double* nsignal_atAT, double* n_pol, double* n_pol_out);
//KD: our rotation matrix
void GetRotation(double theta1, double theta2, double* original_nsignal, double* n_pol, double* n_pol_out);

void GetHitAngle(double*, double*, double*, double*, double&, double&, double&, double&);
void GetHitAngle_ST0(double*, double*, double&, double&, double&, double&);
void GetMirrorHitAngle_ST0(double* posnu2MirrorAT, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component, double& h_component);
void GetHitAngle_ST1(double*, double*, double&, double&, double&, double&, double&, double&, double&, double&);
void GetMirrorHitAngle_ST1(double*, double*, double&, double&, double&, double&, double&, double&, double&, double&);
void GetHitAngle_ST2(double*, double*, double*, double*, double*, double*);
void GetMirrorHitAngle_ST2(double* posnu2MirrorAT, double* n_pol, double* hitangle_e, double* hitangle_h, double* e_component, double* h_component);
void GetHitAngle_LPA(int WhichAntenna, int N_Ant_perST, double* nposnu2AT, double* n_pol, double& hitangle_e, double& hitangle_h, double&, double&); //ST_TYPE=4
//KD modified by adding N_Ant_perST
void GetMirrorHitAngle_LPA(int WhichMirrorAntenna, int N_Ant_perST, double* nposnu2MirrorAT, double* n_pol, double& hitangle_e, double& hitangle_h, double&, double&); //ST_TYPE=4
int Getifreq(double freq);
void GetGZK();//ESS spectrum, 10^16 eV to 10^21.5 eV
void GetGZK7();//energy spectrum starts at 10^17 eV and ends at 10^20 eV
void GetE2();//16-21.5
double PickEnergy(double*, double*); //16-21.5
void GetGZK_fan();//Peter's GZK plot, 17-20
void GetE2_fan();//My E-2 flux, 17-20
double PickEnergy_fan(double*, double*); //17-20

//double GetTauRegen(double , double); //KD 10/25

double  GetGZKIron();

//double Random(double ,double );
void Zero(double*, int);
void Zero(int*, int);
double Distance(double*, double*);
double Square(double);
double Dot(double*, double*);
void Times(double, double*, double*);
double Mag(double*);
void SphAngles(double*, double*);   //KD my spherical angle calculations
double sphangles[2];
double Angle(double*, double*);
void  Cross(double*, double*, double*);
void VectorMinus(double*, double*, double*);
void VectorPlus(double*, double*, double*);
void  nVector(double* v1, double* nv1);
void outVector(double* v, int n);
double Min(double* array, int n);
double Max(double* array, int n);
void GetMaxDis(double pnu, double& Max_distance);
double GetLPM(double);
double Gety();
double fx(double x, double h1, double h2, double nconst, double deltax);

//Fn's for reading in from XML
int SetIntValueXML(tinyxml2::XMLNode * pRoot, int &IntParam, const char* ParamName);
int SetBoolValueXML(tinyxml2::XMLNode * pRoot, int &BoolParam, const char* ParamName);
int SetDoubleValueXML(tinyxml2::XMLNode * pRoot, double &DoubleParam, const char* ParamName);
int SetTextValueXML(tinyxml2::XMLNode * pRoot, const char * &TextParam, const char* ParamName);
int SetElementXML(tinyxml2::XMLNode * pRoot, tinyxml2::XMLElement * &Element, const char* ParamName);
int SetElementXML(tinyxml2::XMLElement * Element, tinyxml2::XMLElement * &SubElement, const char* ParamName);
int ReadStnGeo(const char * infn, int &NAntPerStn, std::vector<AntennaPlacement> &VectAntennas);
int ReadInputXML(const char * infn);


//some const parameters of VmMHz1m
double freq0 = 1150.; //MHz
double f0;
double energy[12], EdNdEdAdt[12]; //KD 10/27 has to be changed for WIDESPECTRUM from 8 to 12
double maxflux;
