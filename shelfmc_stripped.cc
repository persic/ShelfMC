#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <locale>
#include <vector>
#include "TChain.h"
#include "TH1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TStyle.h"
#include "declaration.hh"
#include "TMath.h"
typedef unsigned int uint;
using namespace std;

TStyle* RootStyle();
TStyle* color = RootStyle();

//functions including string type parameters
void GetNuFlavor(string& nuflavor);
void GetCurrent(string& current);
void GetEmHadFrac(string nuflavor, string current, double elast_y, double theta_nu, double L_TauDecay, double L_Ice, double& emfrac, double& hadfrac);
double GetTauRegen(string current, double energy, double theta_nu, double L_TauDecay, double H); //CP 07/15
void CloseTFile(TFile* hfile);
void GetNextNumber(ifstream& in, string& number);
void ReadInput();

//input and output files
//string sgain="/u/cosmic2/dookayka/minna/forLisa/50ohm_horn_gain_down100mhz.txt";
//string sgain="/u/cosmic2/dookayka/minna/forLisa/50ohm_horn_gain.txt";
//string sgain="/u/cosmic/jchanson/NuSim/forLisa/LP_gain_manual.txt";
//string sgain="/pub/jtatar/cshelfmc/forLisa/LP_gain_manual.txt";
// CJR: 2015-07-14 gainsfile moved into ReadGains
//ifstream gainsfile("/pub/jtatar/cshelfmc/forLisa/LP_gain_manual.txt");
ifstream inputfile;
ofstream output;
ofstream outantposall;

TRandom3 Rand3;

//adding flexability to get away from global variables
const int StationType = ST_TYPE;
int main(int argc, char** argv) //MC IceShelf 09/01/2005
{
   string workDir;
   if (argc > 1) {
      workDir = argv[1];
   } else {
      cout << "Expect working directory as input." << endl;
      exit(1);
   }

   string outTag;
   if (argc > 2) {
      outTag = argv[2];
   }
   cout << "Using output tag [" << outTag.c_str() << "]" << endl;

   string input = workDir + "/input.txt";
   inputfile.open(input.c_str());

   string outputfn = workDir + "/output" + outTag + ".txt";
   output.open(outputfn.c_str(), ios::app);

   string outantfn = workDir + "/myantposall_file" + outTag + ".txt";
   outantposall.open(outantfn.c_str(), ios::out);

   time_t raw_start_time = time(NULL);
   struct tm* start_time = localtime(&raw_start_time);

   //cout<<(unsigned)time(NULL);
   //clock_t init_time=clock(); //KD: commented out as it was giving error during compilation
   //output<<asctime(start_time)<<endl;
   //cout<<init_time<<"   "<<asctime(start_time)<<endl;

   cout << asctime(start_time) << endl;

   //Rand3.SetSeed((unsigned)time(NULL));  // CJR: never use time as seed!!

   ReadInput();

   //assigning some globals
   BW = FREQ_HIGH - FREQ_LOW; //MHz
   FREQ_BIN = BW / NFREQ;

   int N_ST_required = 1;
   double Max_distance;

   //ATTEN_UP=(double)atof(argv[1]);
   //REFLECT_RATE=(double)atof(argv[7]);
   //EXPONENT=(double)atof(argv[7]);  //energy exponent
   //JAIME_FACTOR=(double)atof(argv[1]);
   //N_ST_required=(int)atoi(argv[1]);
   //NFIRN=(double)atof(argv[2]);
   //ATGap=(double)atof(argv[8]);
   //SCATTER_WIDTH=atof(argv[1]);

   //ICETHICK=(double)atof(argv[7]);  //KD used rarely, mainly set as a constant
   //ATTEN_UP = 700 - 444*REFLECT_RATE + 200*REFLECT_RATE*REFLECT_RATE; //KD: crude parametrization of curve

   //    NNU = (int)atoi(argv[1]);
   //    N_Ant_perST = (int)atoi(argv[2]); //KD overrides input file number
   //    N_Ant_Trigger = (int)atoi(argv[3]); //KD overrides input file number
   //    NSIGMA = (double)atof(argv[5]); //KD overrides NSIGMA
   //
   //    gainv = (double)atof(argv[6]);
//    gainv = 9.0; //============================JCH June29, 2012
   //gainv = 6.0; // CJR - test
   cout << "hexagonal=" << HEXAGONAL << endl;
   cout << "gainv=" << gainv << endl;
   cout << "signal fluct=" << SIGNAL_FLUCT << endl;
   cout << "tau regen=" << TAUREGENERATION << endl;
   cout << "NNU=" << NNU << endl;
   //
   //    //HRAfactor = (double)atof(argv[8]);
   //    FIRNfactor = (double)atof(argv[7]);
   //    FIRNDEPTH = FIRNDEPTH/FIRNfactor;

   //  ofstream outevents("events.txt",ios::app);
   //  loop over energy
   //  const double EBIN=0.5;
   //  const int LINES=7;
   //  double flux[]={-9.18,-8.5,-7.8,-7.2,-6.775,-6.625,-6.625,-6.75,-7.075,-7.7,-8.8};//16-21
   //  GZK flux, upper blue line in gzkcomp.png in the "events" directory
   //  double flux[]={-7.8,-7.2,-6.775,-6.625,-6.625,-6.75,-7.075,-7.7,-8.8};//17-21
   //  double flux[]={-6.7,-6.7,-6.7,-6.7,-6.7,-6.7,-6.7,-6.7,-6.7,-6.7,-6.7}//E-2 model
   //
   //  for calculating events
   //  double flux[]={-7.8,-7.2,-6.775,-6.625,-6.625,-6.75,-7.075};//1.e17 to 1.e20 eV
   //  double sigma_cm;
   //  double integ=0.;
   //  double constant=9.3623e45;
   //  double dlne=EBIN*log(10.);
   //
   //  for(int ienergy=0;ienergy<LINES;ienergy++)
   //  {
   //  EXPONENT=17.+ienergy*0.5; //loop over energy



   //event information of ST1 for Jiwoo's reconstruction
   typedef struct {
      //Can be filled in before or after ST_TYPE
      Int_t ievt;//which neutrino //neutrino id
      Float_t Energy;//then energy of the neutrino
      Int_t flavor;//the flavor the the neutrino
      Int_t mycurrent;
      Float_t y;//the y factor of the interaction
      Float_t weight;//the weight of the event
      Int_t SumHits;//the sum of all hits including all hits ge than 0.5*NSIGMA
      Int_t N_Ant_perST; //KD: added as a variable
      Int_t NTriggeredST;//number of trigger stations
      Float_t Dir_nu[3];//direction of the neutrino
      Float_t Posi_Int[3];//interaction position

      Float_t entrypt[3];
      Float_t bouncept[3];

      Float_t viewangle_triggered;  //KD: vangle variable
      Float_t viewangle_triggered_mirror; //KD: mirror variable
      Float_t attenfactor;
      Float_t attenfactor_mirror;
      Float_t dis;
      Float_t dis_mirror;
      Float_t vdis;
      Float_t vdis_mirror;
      Float_t cosz;
      Float_t cosz_mirror;
      Float_t coszp;
      Float_t coszp_mirror;

      Float_t theta_signal;
      Float_t theta_signal_atAT;
      Float_t theta_signal_mirror;
      Float_t theta_signal_atAT_mirror;
      Float_t phi_nposnu2ST;
      Float_t phi_nposnu2MirrorST;
      Float_t theta_nposnu2ST;
      Float_t theta_nposnu2MirrorST;

      Float_t my_time_AT2ST[100];
      Float_t my_time_AT2ST_mirror[100];
      Float_t nsignal_i_atAT[100];
      Float_t nsignal_j_atAT[100];
      Float_t nsignal_k_atAT[100];
      Float_t nposnu2AT_i[100];
      Float_t nposnu2AT_j[100];
      Float_t nposnu2AT_k[100];
      Float_t nsignal_mirror_i_atAT[100];
      Float_t nsignal_mirror_j_atAT[100];
      Float_t nsignal_mirror_k_atAT[100];
      Float_t nposnu2MirrorAT_i[100];
      Float_t nposnu2MirrorAT_j[100];
      Float_t nposnu2MirrorAT_k[100];

      Float_t my_ST_Posi[2][3];
      Float_t my_ST_abs_time;
      Float_t my_ST_abs_time_mirror;
      Float_t my_d_posnu2ST;
      Float_t my_d_posnu2MirrorST;
      Float_t my_nsignal_atST[2][3];

      Float_t Bfield_D_firn[3];
      Float_t Bfield_D_original[3];
      Float_t Polarization_D[3];
      Float_t Polarization_R[3];
      Float_t polariz_theta_D;
      Float_t polariz_phi_D;
      Float_t polariz_theta_R;
      Float_t polariz_phi_R;

      Float_t Fresnel_Pol_D[3];
      Float_t Fresnel_Pol_D_phi;
      Float_t Fresnel_Pol_D_theta;
      Float_t Fresnel_Pol_R[3];
      Float_t Fresnel_Pol_R_phi;
      Float_t Fresnel_Pol_R_theta;

      //at vertex now
      Float_t original_nsignal[3];
      Float_t original_nsignal_mirror[3];
      Float_t Polarization_D_original[3];
      Float_t Polarization_R_original[3];
      Float_t polariz_theta_D_original;
      Float_t polariz_phi_D_original;
      Float_t polariz_theta_R_original;
      Float_t polariz_phi_R_original;

      //these here are for straight line from interaction point to station, useful for no firn case. CAUTION otherwise.
      Float_t Polariz_reflect[3];
      Float_t Polariz_direct[3];
      Float_t Bfield_direct[3];
      Float_t Bfield_reflect[3];
      Float_t Polariz_direct_phi;
      Float_t Polariz_direct_theta;
      Float_t Polariz_reflect_phi;
      Float_t Polariz_reflect_theta;

      Float_t theta_nu;          //KD: hh1 histogram
      Float_t costhetanu;        //KD: costhetanu histogram
      Float_t depth;          //KD: depth histogram
      Float_t phi_nu;

      Int_t sum_triggeredST;
      Int_t sum_triggeredST_mirror;

      Int_t iAT[100];   //KD labels AT within a station
      Int_t iAT_mirror[100];
      Int_t iAT_index_mymax; //corresponds to voltage
      Int_t iAT_mirror_index_mymax;  //corresponds to voltage
      Int_t iAT_index_mymin; //corresponds to voltage
      Int_t iAT_mirror_index_mymin;  //corresponds to voltage
      Int_t iAT_index_timemax; //corresponds to time
      Int_t iAT_index_timemin; //corresponds to time
      Int_t iAT_mirror_index_timemax; //corresponds to time
      Int_t iAT_mirror_index_timemin; //corresponds to time

      Float_t volt_LPA[100];
      Float_t volt_LPA_preNoise[100];
      Float_t hitangle_e_LPA[100];
      Float_t hitangle_h_LPA[100];
      Float_t theta_my_signalAT[100];
      Float_t phi_my_signalAT[100];
      Float_t hitangle_e_LPA_mymax;
      Float_t hitangle_h_LPA_mymax;
      Float_t totaltime[100];
      Float_t totaltime_mymax;
      Float_t totaltime_mymin;

      //----------------------
      //COMMENTED OUT ON 5/12/11
      Float_t term_LPA_my[100][95];
      Float_t heff_my[95];
      Float_t vmmhz_my[100][95];
      Float_t vmmhz_mirror_my[100][95];
      Float_t term_LPA_mirror_my[100][95];
      Float_t Efield_atposnu_my[100][95];
      Float_t Efield_atposnu_mirror_my[100][95];
      //------------------------

      Float_t freq_my[95];    //re-include in b1 for angular reconstruction

      Float_t sum_vmmhz1m_unattened_max;
      Float_t Efield_atposnu_mymax;
      Float_t volte_allAT_mymax;

      //JCH: Added this one little variable
      Float_t volte_allAT_mymin;

      Float_t volt_LPA_mirror[100];
      Float_t volt_LPA_mirror_preNoise[100];
      Float_t hitangle_e_LPA_mirror[100];
      Float_t hitangle_h_LPA_mirror[100];
      Float_t theta_my_signalAT_mirror[100];
      Float_t phi_my_signalAT_mirror[100];
      Float_t hitangle_e_LPA_mirror_mymax;
      Float_t hitangle_h_LPA_mirror_mymax;
      Float_t totaltime_mirror[100];
      Float_t totaltime_mirror_mymax;
      Float_t totaltime_mirror_mymin;

      Float_t sum_vmmhz1m_unattened_mirror_max;
      Float_t Efield_atposnu_mirror_mymax;
      Float_t volte_allAT_mirror_mymax;

      Float_t elpm;
      Float_t vmmhz1m_max;
      Float_t changle;
      Float_t eshower_em;
      Float_t eshower_had;
      Float_t e_component_LPA[100];
      Float_t h_component_LPA[100];
      Float_t e_component_LPA_mirror[100];
      Float_t h_component_LPA_mirror[100];

   } Event_Info;

   Event_Info b1;


   //COMMENTED OUT ON 5/11/11 for lower output file
   //==============================================================
   //creating a tree b7 that contains Energy tracking variables
   struct {
      Int_t N_Ant_perST; //KD: added as a variable
      Float_t Evmmhz1m_max;
      Float_t Edeltheta_em_max;
      Float_t Edeltheta_had_max;
      Float_t En1;
      Float_t Eshowerlength;
      Float_t Evmmhz[100][95];
      Float_t Evmmhz_mirror[100][95];
      Float_t Esum_vmmhz[100];
      Float_t Esum_vmmhz_mirror[100];
      Float_t Evmmhz1m_unattened[95];
      Float_t Evmmhz1m_unattened_mirror[95];
      Float_t Esumvmmhz1m_unattened;
      Float_t Esumvmmhz1m_unattened_mirror;
      Float_t Edeltheta_em[95];
      Float_t Edeltheta_had[95];
      Float_t Eemfrac;
      Float_t Ehadfrac;
      Float_t Evmmhz_taper[100][95];
      Float_t Evmmhz_taper_mirror[100][95];
      Float_t Esum_vmmhz_taper[100];
      Float_t Esum_vmmhz_taper_mirror[100];
      Float_t Evmmhz_em[100][95];
      Float_t Evmmhz_em_mirror[100][95];
      Float_t Evmmhz_had_mirror[100][95];
      Float_t Evmmhz_had[100][95];

      //07/17/11 shifting some variables from b1 into here
      Float_t term_LPA_my[100][95];
      Float_t heff_my[95];
      Float_t vmmhz_my[100][95];
      Float_t freq_my[95];  //kept a copy in b1.
      Float_t vmmhz_mirror_my[100][95];
      Float_t term_LPA_mirror_my[100][95];
      Float_t Efield_atposnu_my[100][95];
      Float_t Efield_atposnu_mirror_my[100][95];

   } b7;



   //======================================================================================+
   //======================== START OF BUILDING ROOT FILES ================================+
   //======================================================================================+
   string file2 = workDir + "/ShelfMCTrees" + outTag + ".root";
   TFile* hfile = new TFile(file2.c_str(), "RECREATE", "iceshelf"); //KD: shifted above.


   //------------------------------------------------------------------+
   TTree* tree1 = new TTree("PAM", "a tree of triggered events"); //    |
   //------------------------------------------------------------------+

   tree1->Branch("ievt", &b1.ievt, "ievt/I");
   tree1->Branch("Energy", &b1.Energy, "Energy/F");
   tree1->Branch("flavor", &b1.flavor, "flavor/I");
   tree1->Branch("mycurrent", &b1.mycurrent, "mycurrent/I");
   tree1->Branch("y", &b1.y, "y/F");
   tree1->Branch("weight", &b1.weight, "weight/F");
   tree1->Branch("SumHits", &b1.SumHits, "SumHits/I");
   //tree1->Branch("N_allchannels",&b1.N_allchannels,"N_allchannels/I");
   tree1->Branch("N_Ant_perST", &b1.N_Ant_perST, "N_Ant_perST/I");
   tree1->Branch("NTriggeredST", &b1.NTriggeredST, "NTriggeredST/I");
   tree1->Branch("Dir_nu", b1.Dir_nu, "Dir_nu[3]/F");
   tree1->Branch("Posi_Int", b1.Posi_Int, "Posi_Int[3]/F");

   tree1->Branch("entrypt", b1.entrypt, "entrypt[3]/F");
   tree1->Branch("bouncept", b1.bouncept, "bouncept[3]/F");

   //tree1->Branch("Trig_Type",b1.Trig_Type,"Trig_Type[SumHits]/I");
   //tree1->Branch("Abs_Time",b1.Abs_Time,"Abs_Time[SumHits]/F");
   //tree1->Branch("Timedelay",b1.Timedelay,"Timedelay[SumHits]/D");
   //tree1->Branch("Path_inice",b1.Path_inice,"Path_inice[SumHits]/D");
   //tree1->Branch("Path_infirn",b1.Path_infirn,"Path_infirn[SumHits]/D");
   //tree1->Branch("Amp",b1.Amp,"Amp[SumHits]/F");
   //tree1->Branch("E_atposnu",b1.E_atposnu,"E_atposnu[SumHits]/D");
   //tree1->Branch("Emax_atposnu",b1.Emax_atposnu,"Emax_atposnu[SumHits]/D");
   //tree1->Branch("N_channels_fired",b1.N_channels_fired,"N_channels_fired[SumHits]/I");
   //tree1->Branch("Volts",b1.Volts,"Volts[N_allchannels]/D");
   //tree1->Branch("ST_Posi_x",b1.ST_Posi_x,"ST_Posi_x[SumHits]/F");
   //tree1->Branch("ST_Posi_y",b1.ST_Posi_y,"ST_Posi_y[SumHits]/F");
   //tree1->Branch("ST_Posi_z",b1.ST_Posi_z,"ST_Posi_z[SumHits]/F");
   //tree1->Branch("ST_iRow",b1.ST_iRow,"ST_iRow[SumHits]/I");
   //tree1->Branch("ST_iCol",b1.ST_iCol,"ST_iCol[SumHits]/I");
   //tree1->Branch("iSIGMA",b1.iSIGMA,"iSIGMA[SumHits]/I");
   //tree1->Branch("Dir_Shower",b1.Dir_Shower,"Dir_Shower[SumHits]/F");

   //KD: my adding new variables here
   tree1->Branch("viewangle_triggered", &b1.viewangle_triggered, "viewangle_triggered/F");
   tree1->Branch("viewangle_triggered_mirror", &b1.viewangle_triggered_mirror, "viewangle_triggered_mirror/F");
   tree1->Branch("attenfactor", &b1.attenfactor, "attenfactor/F");
   tree1->Branch("attenfactor_mirror", &b1.attenfactor_mirror, "attenfactor_mirror/F");
   tree1->Branch("dis", &b1.dis, "dis/F");
   tree1->Branch("dis_mirror", &b1.dis_mirror, "dis_mirror/F");
   tree1->Branch("vdis", &b1.vdis, "vdis/F");
   tree1->Branch("vdis_mirror", &b1.vdis_mirror, "vdis_mirror/F");
   tree1->Branch("cosz", &b1.cosz, "cosz/F");
   tree1->Branch("cosz_mirror", &b1.cosz_mirror, "cosz_mirror/F");
   tree1->Branch("coszp", &b1.coszp, "coszp/F");
   tree1->Branch("coszp_mirror", &b1.coszp_mirror, "coszp_mirror/F");

   tree1->Branch("theta_signal", &b1.theta_signal, "theta_signal/F");
   tree1->Branch("theta_signal_atAT", &b1.theta_signal_atAT, "theta_signal_atAT/F");
   tree1->Branch("theta_signal_mirror", &b1.theta_signal_mirror, "theta_signal_mirror/F");
   tree1->Branch("theta_signal_atAT_mirror", &b1.theta_signal_atAT_mirror, "theta_signal_atAT_mirror/F");
   tree1->Branch("phi_nposnu2ST", &b1.phi_nposnu2ST, "phi_nposnu2ST(in deg)/F");
   tree1->Branch("phi_nposnu2MirrorST", &b1.phi_nposnu2MirrorST, "phi_nposnu2MirrorST(in deg)/F");
   tree1->Branch("theta_nposnu2ST", &b1.theta_nposnu2ST, "theta_nposnu2ST(in deg)(CAUTION)/F");
   tree1->Branch("theta_nposnu2MirrorST", &b1.theta_nposnu2MirrorST, "theta_nposnu2MirrorST(CAUTION)(in deg)/F");

   //tree1->Branch("nposnu2ST_i", &b1.nposnu2ST_i, "nposnu2ST_i/F");
   //tree1->Branch("nposnu2ST_j", &b1.nposnu2ST_j, "nposnu2ST_j/F");
   //tree1->Branch("nposnu2ST_k", &b1.nposnu2ST_k, "nposnu2ST_k/F");

   tree1->Branch("my_time_AT2ST", &b1.my_time_AT2ST, "my_time_AT2ST[N_Ant_perST]/F");
   tree1->Branch("my_time_AT2ST_mirror", &b1.my_time_AT2ST_mirror, "my_time_AT2ST_mirror[N_Ant_perST]/F");
   tree1->Branch("nsignal_i_atAT", &b1.nsignal_i_atAT, "nsignal_i_atAT[N_Ant_perST]/F");
   tree1->Branch("nsignal_j_atAT", &b1.nsignal_j_atAT, "nsignal_j_atAT[N_Ant_perST]/F");
   tree1->Branch("nsignal_k_atAT", &b1.nsignal_k_atAT, "nsignal_k_atAT[N_Ant_perST]/F");
   tree1->Branch("nposnu2AT_i", &b1.nposnu2AT_i, "nposnu2AT_i[N_Ant_perST](CAUTION)/F");
   tree1->Branch("nposnu2AT_j", &b1.nposnu2AT_j, "nposnu2AT_j[N_Ant_perST](CAUTION)/F");
   tree1->Branch("nposnu2AT_k", &b1.nposnu2AT_k, "nposnu2AT_k[N_Ant_perST](CAUTION)/F");
   tree1->Branch("nsignal_mirror_i_atAT", &b1.nsignal_mirror_i_atAT, "nsignal_mirror_i_atAT[N_Ant_perST]/F");
   tree1->Branch("nsignal_mirror_j_atAT", &b1.nsignal_mirror_j_atAT, "nsignal_mirror_j_atAT[N_Ant_perST]/F");
   tree1->Branch("nsignal_mirror_k_atAT", &b1.nsignal_mirror_k_atAT, "nsignal_mirror_k_atAT[N_Ant_perST]/F");
   tree1->Branch("nposnu2MirrorAT_i", &b1.nposnu2MirrorAT_i, "nposnu2MirrorAT_i[N_Ant_perST]/F");
   tree1->Branch("nposnu2MirrorAT_j", &b1.nposnu2MirrorAT_j, "nposnu2MirrorAT_j[N_Ant_perST]/F");
   tree1->Branch("nposnu2MirrorAT_k", &b1.nposnu2MirrorAT_k, "nposnu2MirrorAT_k[N_Ant_perST]/F");

   tree1->Branch("my_nsignal_atST", &b1.my_nsignal_atST, "my_nsignal_atST[2][3](both DIR & REF)/F");
   tree1->Branch("my_ST_Posi", &b1.my_ST_Posi, "my_ST_Posi[2][3]/F");
   tree1->Branch("my_ST_abs_time", &b1.my_ST_abs_time, "my_ST_abs_time/F");
   tree1->Branch("my_ST_abs_time_mirror", &b1.my_ST_abs_time_mirror, "my_ST_abs_time_mirror/F");
   tree1->Branch("my_d_posnu2ST", &b1.my_d_posnu2ST, "my_d_posnu2ST/F");
   tree1->Branch("my_d_posnu2MirrorST", &b1.my_d_posnu2MirrorST, "my_d_posnu2MirrorST/F");

   tree1->Branch("Polarization_D", &b1.Polarization_D, "Polarization_D[3](at AT)/F");
   tree1->Branch("Bfield_D_firn", &b1.Bfield_D_firn, "Bfield_D_firn[3](at AT)/F");

   tree1->Branch("Polarization_R", &b1.Polarization_R, "Polarization_R[3](at AT)/F");
   tree1->Branch("polariz_theta_D", &b1.polariz_theta_D, "polariz_theta_D(in deg)(at AT)/F");
   tree1->Branch("polariz_phi_D", &b1.polariz_phi_D, "polariz_phi_D(in deg)(at AT)/F");
   tree1->Branch("polariz_theta_R", &b1.polariz_theta_R, "polariz_theta_R(in deg)(at AT)/F");
   tree1->Branch("polariz_phi_R", &b1.polariz_phi_R, "polariz_phi_R(in deg)(at AT)/F");

   tree1->Branch("Fresnel_Pol_D", &b1.Fresnel_Pol_D, "Fresnel_Pol_D[3]/F");
   tree1->Branch("Fresnel_Pol_D_phi", &b1.Fresnel_Pol_D_phi, "Fresnel_Pol_D_phi/F");
   tree1->Branch("Fresnel_Pol_D_theta", &b1.Fresnel_Pol_D_theta, "Fresnel_Pol_D_theta/F");
   tree1->Branch("Fresnel_Pol_R", &b1.Fresnel_Pol_R, "Fresnel_Pol_R[3]/F");
   tree1->Branch("Fresnel_Pol_R_phi", &b1.Fresnel_Pol_R_phi, "Fresnel_Pol_R_phi/F");
   tree1->Branch("Fresnel_Pol_R_theta", &b1.Fresnel_Pol_R_theta, "Fresnel_Pol_R_theta/F");


   tree1->Branch("original_nsignal", &b1.original_nsignal, "original_nsignal[3]/F");
   tree1->Branch("original_nsignal_mirror", &b1.original_nsignal_mirror, "original_nsignal_mirror[3]/F");
   tree1->Branch("Polarization_D_original", &b1.Polarization_D_original, "Polarization_D_original[3](at vertex)/F");
   tree1->Branch("Bfield_D_original", &b1.Bfield_D_original, "Bfield_D_original[3](at vertex)/F");
   tree1->Branch("Polarization_R_original", &b1.Polarization_R_original, "Polarization_R_original[3](at vertex)/F");
   tree1->Branch("polariz_theta_D_original", &b1.polariz_theta_D_original, "polariz_theta_D_original(in deg)(at vertex)/F");
   tree1->Branch("polariz_phi_D_original", &b1.polariz_phi_D_original, "polariz_phi_D_original(in deg)(at vertex)/F");
   tree1->Branch("polariz_theta_R_original", &b1.polariz_theta_R_original, "polariz_theta_R_original(in deg)(at vertex)/F");
   tree1->Branch("polariz_phi_R_original", &b1.polariz_phi_R_original, "polariz_phi_R_original(in deg)(at vertex)/F");



   tree1->Branch("Polariz_direct", &b1.Polariz_direct, "Polariz_direct[3](original)(CAUTION)/F");
   tree1->Branch("Polariz_reflect", &b1.Polariz_reflect, "Polariz_reflect[3](original)(CAUTION)/F");
   tree1->Branch("Bfield_direct", &b1.Bfield_direct, "Bfield_direct[3](original)(CAUTION)/F");
   tree1->Branch("Bfield_reflect", &b1.Bfield_reflect, "Bfield_reflect[3](original)(CAUTION)/F");

   tree1->Branch("Polariz_direct_phi", &b1.Polariz_direct_phi, "Polariz_direct_phi(original, in deg)(CAUTION)/F");
   tree1->Branch("Polariz_direct_theta", &b1.Polariz_direct_theta, "Polariz_direct_theta(original, in deg)(CAUTION)/F");
   tree1->Branch("Polariz_reflect_phi", &b1.Polariz_reflect_phi, "Polariz_reflect_phi(original, in deg)(CAUTION)/F");
   tree1->Branch("Polariz_reflect_theta", &b1.Polariz_reflect_theta, "Polariz_reflect_theta(original, in deg)(CAUTION)/F");

   tree1->Branch("theta_nu", &b1.theta_nu, "theta_nu (in degrees)/F");
   tree1->Branch("costhetanu", &b1.costhetanu, "costhetanu/F");
   tree1->Branch("depth", &b1.depth, "depth/F");
   tree1->Branch("phi_nu", &b1.phi_nu, "phi_nu (in degrees)/F");

   tree1->Branch("sum_triggeredST", &b1.sum_triggeredST, "sum_triggeredST/I");
   tree1->Branch("sum_triggeredST_mirror", &b1.sum_triggeredST_mirror, "sum_triggeredST_mirror/I");

   tree1->Branch("iAT", &b1.iAT, "iAT[N_Ant_perST](NOTE +1)/I");
   tree1->Branch("iAT_mirror", &b1.iAT_mirror, "iAT_mirror[N_Ant_perST](NOTE +1)/I");
   tree1->Branch("iAT_index_mymax", &b1.iAT_index_mymax, "iAT_index_mymax/I");
   tree1->Branch("iAT_mirror_index_mymax", &b1.iAT_mirror_index_mymax, "iAT_mirror_index_mymax/I");
   tree1->Branch("iAT_index_mymin", &b1.iAT_index_mymin, "iAT_index_mymin/I");
   tree1->Branch("iAT_mirror_index_mymin", &b1.iAT_mirror_index_mymin, "iAT_mirror_index_mymin/I");
   tree1->Branch("iAT_index_timemax", &b1.iAT_index_timemax, "iAT_index_timemax/I");
   tree1->Branch("iAT_mirror_index_timemax", &b1.iAT_mirror_index_timemax, "iAT_mirror_index_timemax/I");
   tree1->Branch("iAT_index_timemin", &b1.iAT_index_timemin, "iAT_index_timemin/I");
   tree1->Branch("iAT_mirror_index_timemin", &b1.iAT_mirror_index_timemin, "iAT_mirror_index_timemin/I");
   //tree1->Branch("N_Ant_perST", &b1.N_Ant_perST, "N_Ant_perST/I");
   //tree1->Branch("kviewangle", &b1.kviewangle, "kviewangle[N_Ant_perST]/F");
   //tree1->Branch("volt_max", &b1.volt_max, "volt_max[N_Ant_perST]/F");
   tree1->Branch("volt_LPA", &b1.volt_LPA, "volt_LPA[N_Ant_perST]/F");
   tree1->Branch("volt_LPA_preNoise", &b1.volt_LPA_preNoise, "volt_LPA_preNoise[N_Ant_perST]/F");
   //tree1->Branch("sum_vmmhz", &b1.sum_vmmhz, "sum_vmmhz[N_Ant_perST]/F");
   //tree1->Branch("sum_vmmhz1m_unattened", &b1.sum_vmmhz1m_unattened, "sum_vmmhz1m_unattened[N_Ant_perST]/F");
   tree1->Branch("hitangle_e_LPA", &b1.hitangle_e_LPA, "hitangle_e_LPA[N_Ant_perST]/F");
   tree1->Branch("hitangle_h_LPA", &b1.hitangle_h_LPA, "hitangle_h_LPA[N_Ant_perST]/F");
   tree1->Branch("theta_my_signalAT", &b1.theta_my_signalAT, "theta_my_signalAT[N_Ant_perST]/F");
   tree1->Branch("phi_my_signalAT", &b1.phi_my_signalAT, "phi_my_signalAT[N_Ant_perST]/F");
   tree1->Branch("hitangle_e_LPA_mymax", &b1.hitangle_e_LPA_mymax, "hitangle_e_LPA_mymax/F");
   tree1->Branch("hitangle_h_LPA_mymax", &b1.hitangle_h_LPA_mymax, "hitangle_h_LPA_mymax/F");
   tree1->Branch("totaltime", &b1.totaltime, "totaltime[N_Ant_perST]/F");
   tree1->Branch("totaltime_mymin", &b1.totaltime_mymin, "totaltime_mymin/F");
   tree1->Branch("totaltime_mymax", &b1.totaltime_mymax, "totaltime_mymax/F");
   tree1->Branch("sum_vmmhz1m_unattened_max", &b1.sum_vmmhz1m_unattened_max, "sum_vmmhz1m_unattened_max/F");
   tree1->Branch("Efield_atposnu_mymax", &b1.Efield_atposnu_mymax, "Efield_atposnu_mymax/F");
   tree1->Branch("volte_allAT_mymax", &b1.volte_allAT_mymax, "volte_allAT_mymax/F");
   tree1->Branch("volte_allAT_mymin", &b1.volte_allAT_mymin, "volte_allAT_mymin/F");

   tree1->Branch("hitangle_e_LPA_mirror", &b1.hitangle_e_LPA_mirror, "hitangle_e_LPA_mirror[N_Ant_perST]/F");
   tree1->Branch("hitangle_h_LPA_mirror", &b1.hitangle_h_LPA_mirror, "hitangle_h_LPA_mirror[N_Ant_perST]/F");
   tree1->Branch("theta_my_signalAT_mirror", &b1.theta_my_signalAT_mirror, "theta_my_signalAT_mirror[N_Ant_perST]/F");
   tree1->Branch("phi_my_signalAT_mirror", &b1.phi_my_signalAT_mirror, "phi_my_signalAT_mirror[N_Ant_perST]/F");
   tree1->Branch("hitangle_e_LPA_mirror_mymax", &b1.hitangle_e_LPA_mirror_mymax, "hitangle_e_LPA_mirror_mymax/F");
   tree1->Branch("hitangle_h_LPA_mirror_mymax", &b1.hitangle_h_LPA_mirror_mymax, "hitangle_h_LPA_mirror_mymax/F");
   tree1->Branch("totaltime_mirror", &b1.totaltime_mirror, "totaltime_mirror[N_Ant_perST]/F");
   tree1->Branch("totaltime_mirror_mymin", &b1.totaltime_mirror_mymin, "totaltime_mirror_mymin/F");
   tree1->Branch("totaltime_mirror_mymax", &b1.totaltime_mirror_mymax, "totaltime_mirror_mymax/F");

   //-----------------------
   // COMMENTED OUT ON 5/12/11 //SHIFTED TO B7 ON 07/17/11
   tree1->Branch("term_LPA_my", &b1.term_LPA_my, "term_LPA_my[N_Ant_perST][95]/F");
   tree1->Branch("heff_my", &b1.heff_my, "heff_my[95]/F");
   tree1->Branch("vmmhz_my", &b1.vmmhz_my, "vmmhz_my[N_Ant_perST][95]/F");
   tree1->Branch("vmmhz_mirror_my", &b1.vmmhz_mirror_my, "vmmhz_mirror_my[N_Ant_perST][95]/F");
   tree1->Branch("term_LPA_mirror_my", &b1.term_LPA_mirror_my, "term_LPA_mirror_my[N_Ant_perST][95]/F");
   tree1->Branch("Efield_atposnu_my", &b1.Efield_atposnu_my, "Efield_atposnu_my[N_Ant_perST][95]/F");
   tree1->Branch("Efield_atposnu_mirror_my", &b1.Efield_atposnu_mirror_my, "Efield_atposnu_mirror_my[N_Ant_perST][95]/F");
   //-----------------------
   tree1->Branch("freq_my", &b1.freq_my, "freq_my[95]/F"); //re-included as in b7

   //tree1->Branch("kviewangle_mirror", &b1.kviewangle_mirror, "kviewangle_mirror[N_Ant_perST]/F");
   //tree1->Branch("volt_mirror_max", &b1.volt_mirror_max, "volt_mirror_max[N_Ant_perST]/F");
   tree1->Branch("volt_LPA_mirror", &b1.volt_LPA_mirror, "volt_LPA_mirror[N_Ant_perST]/F");
   tree1->Branch("volt_LPA_mirror_preNoise", &b1.volt_LPA_mirror_preNoise, "volt_LPA_mirror_preNoise[N_Ant_perST]/F");
   //tree1->Branch("sum_vmmhz_mirror", &b1.sum_vmmhz_mirror, "sum_vmmhz_mirror[N_Ant_perST]/F");
   //tree1->Branch("sum_vmmhz1m_unattened_mirror", &b1.sum_vmmhz1m_unattened_mirror, "sum_vmmhz1m_unattened_mirror[N_Ant_perST]/F");
   tree1->Branch("sum_vmmhz1m_unattened_mirror_max", &b1.sum_vmmhz1m_unattened_max, "sum_vmmhz1m_unattened_mirror_max/F");
   tree1->Branch("Efield_atposnu_mirror_mymax", &b1.Efield_atposnu_mirror_mymax, "Efield_atposnu_mirror_mymax/F");
   tree1->Branch("volte_allAT_mirror_mymax", &b1.volte_allAT_mirror_mymax, "volte_allAT_mirror_mymax/F");

   tree1->Branch("elpm", &b1.elpm, "elpm/F");
   tree1->Branch("vmmhz1m_max", &b1.vmmhz1m_max, "vmmhz1m_max/F");
   tree1->Branch("changle", &b1.changle, "changle/F");
   tree1->Branch("eshower_em", &b1.eshower_em, "eshower_em/F");
   tree1->Branch("eshower_had", &b1.eshower_had, "eshower_had/F");

   tree1->Branch("e_component_LPA", &b1.e_component_LPA, "e_component_LPA[N_Ant_perST]/F");
   tree1->Branch("h_component_LPA", &b1.h_component_LPA, "h_component_LPA[N_Ant_perST]/F");
   tree1->Branch("e_component_LPA_mirror", &b1.e_component_LPA_mirror, "e_component_LPA_mirror[N_Ant_perST]/F");
   tree1->Branch("h_component_LPA_mirror", &b1.h_component_LPA_mirror, "h_component_LPA_mirror[N_Ant_perST]/F");


//---------------------------------------------------------------------------------------
// below is an attempt at events info for multi-stations, from original code for ST3
   struct {
      Int_t ievt;
      Float_t Energy;
      Int_t flavor;
      Float_t y;
      Float_t weight;
      Int_t NTriggeredST;
      Float_t Theta_nu;
      Int_t N_ATs;//total number of triggered antennas of each triggered event
      Float_t Dir_nu[3];//direction of neutrino
      Float_t Posi_Int[3];//interaction position

      Int_t NAT_atST[500];//number of the triggered antennas at each triggered station
      Int_t Trig_Type[500];//station triggered by dir, ref or both
      Int_t ST_iRow[500];//row number of the triggered station
      Int_t ST_iCol[500];//col number of the triggered station

      //struct{
      Int_t iAT[4000];//which antenna is being triggered
      //Int_t iAT_mirror[4000];//which antenna is being triggered

      Int_t sum_triggeredST;
      Int_t sum_triggeredST_mirror;
      Int_t nrows;
      Int_t ncols;

      Double_t ST_Posi_x[31][31];//x component
      Double_t ST_Posi_y[31][31];//y component
      Double_t ST_Posi_z[31][31];//z component

      Double_t path_inice[4000];
      Double_t path_infirn[4000];
      Double_t totaltime[4000];//the time to get to each triggered antenna
      Double_t Efield_atposnu[4000];
      Double_t volte_allAT[4000];
      Double_t volth_allAT[4000];
      Float_t AT_Posi_x[4000];//x component of the triggered antenna
      Float_t AT_Posi_y[4000];//y component
      Float_t AT_Posi_z[4000];//z component
      Float_t Dir_shower_x[4000];//x component of the signal reaching the antenna
      Float_t Dir_shower_y[4000];//y component of the signal reaching the antenna
      Float_t Dir_shower_z[4000];//z component of the signal reaching the antenna
      Double_t theta_shower_atposnu[4000];
      // } a1[500];
   } b3;

//------------------------------------------------------------------+
   TTree* tree3 = new TTree("PAMst", "a tree for multi-stations"); //    |
//------------------------------------------------------------------+
   //TTree *tree3=new TTree("ST_TYPE3","a tree for ST_TYPE 3");

   tree3->Branch("ievt", &b3.ievt, "ievt/I");
   tree3->Branch("Energy", &b3.Energy, "Energy/F");
   tree3->Branch("flavor", &b3.flavor, "flavor/I");
   tree3->Branch("y", &b3.y, "y/F");
   tree3->Branch("weight", &b3.weight, "weight/F");
   tree3->Branch("NTriggeredST", &b3.NTriggeredST, "NTriggeredST/I");
   tree3->Branch("Theta_nu", &b3.Theta_nu, "Theta_nu/F");
   tree3->Branch("N_ATs", &b3.N_ATs, "N_ATs/I");
   tree3->Branch("Dir_nu", &b3.Dir_nu, "Dir_nu[3]/F");
   tree3->Branch("Posi_Int", &b3.Posi_Int, "Posi_Int[3]/F");
   tree3->Branch("NAT_atST", &b3.NAT_atST, "NAT_atST[NTriggeredST]/I");
   tree3->Branch("Trig_Type", &b3.Trig_Type, "Trig_Type[NTriggeredST]/I");
   tree3->Branch("ST_iRow", &b3.ST_iRow, "ST_iRow[NTriggeredST]/I");
   tree3->Branch("ST_iCol", &b3.ST_iCol, "ST_iCol[NTriggeredST]/I");

   tree3->Branch("iAT", &b3.iAT, "iAT[N_ATs]/I");
   //tree3->Branch("iAT_mirror",&b3.iAT_mirror,"iAT_mirror[N_ATs]/I");
   tree3->Branch("sum_triggeredST", &b3.sum_triggeredST, "sum_triggeredST/I");
   tree3->Branch("sum_triggeredST_mirror", &b3.sum_triggeredST_mirror, "sum_triggeredST_mirror/I");
   tree3->Branch("nrows", &b3.nrows, "nrows/I");
   tree3->Branch("ncols", &b3.ncols, "ncols/I");
   tree3->Branch("ST_Posi_x", &b3.ST_Posi_x, "ST_Posi_x[31][31]/D");
   tree3->Branch("ST_Posi_y", &b3.ST_Posi_y, "ST_Posi_y[31][31]/D");
   tree3->Branch("ST_Posi_z", &b3.ST_Posi_z, "ST_Posi_z[31][31]/D");

   tree3->Branch("path_inice", &b3.path_inice, "path_inice[N_ATs]/D");
   tree3->Branch("path_infirn", &b3.path_infirn, "path_infirn[N_ATs]/D");
   tree3->Branch("totaltime", &b3.totaltime, "totaltime[N_ATs]/D");
   tree3->Branch("Efield_atposnu", &b3.Efield_atposnu, "Efield_atposnu[N_ATs]/D");
   tree3->Branch("volte_allAT", &b3.volte_allAT, "volte_allAT[N_ATs]/D");
   tree3->Branch("volth_allAT", &b3.volth_allAT, "volth_allAT[N_ATs]/D");
   tree3->Branch("AT_Posi_x", &b3.AT_Posi_x, "AT_Posi_x[N_ATs]/F");
   tree3->Branch("AT_Posi_y", &b3.AT_Posi_y, "AT_Posi_y[N_ATs]/F");
   tree3->Branch("AT_Posi_z", &b3.AT_Posi_z, "AT_Posi_z[N_ATs]/F");
   tree3->Branch("Dir_shower_x", &b3.Dir_shower_x, "Dir_shower_x[N_ATs]/F");
   tree3->Branch("Dir_shower_y", &b3.Dir_shower_y, "Dir_shower_y[N_ATs]/F");
   tree3->Branch("Dir_shower_z", &b3.Dir_shower_z, "Dir_shower_z[N_ATs]/F");
   tree3->Branch("theta_shower_atposnu", &b3.theta_shower_atposnu, "theta_shower_atposnu[N_ATs]/D");


   b3.nrows = NROWS;
   b3.ncols = NCOLS;



//------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------

//creating a tree that specifies the simulation conditions
   struct {
      Int_t nnu;
      Double_t atgap;
      Int_t nrows;
      Int_t ncols;

      Bool_t signal_fluctuation;
      Bool_t lpm_effect;
      Bool_t tau_regeneration;
      Bool_t planewaves;
      Bool_t depthdependentfirn;
      Bool_t firn_presence;
      Bool_t spectrum;
      Bool_t widespectrum;
      Bool_t gzk_flux;
      Bool_t seckel_effect;
      Int_t atten_freq;
      Int_t atten_exp;
      Int_t const_attenlength;
      Int_t shadowing_dir_ref;
   } b5;

//------------------------------------------------------------------+
   TTree* tree5 = new TTree("sim_parameters", "a tree of parameters boolean"); //
//------------------------------------------------------------------+
   tree5->Branch("nnu", &b5.nnu, "nnu/I");
   tree5->Branch("atgap", &b5.atgap, "atgap/D");
   tree5->Branch("nrows", &b5.nrows, "nrows/I");
   tree5->Branch("ncols", &b5.ncols, "ncols/I");
//tree5->Branch("ST_Posi_x", &b5.ST_Posi_x,"ST_Posi_x[nrows][ncols]/D");
//tree5->Branch("ST_Posi_y", &b5.ST_Posi_y,"ST_Posi_y[nrows][ncols]/D");
//tree5->Branch("ST_Posi_z", &b5.ST_Posi_z,"ST_Posi_z[nrows][ncols]/D");

   tree5->Branch("signal_fluctuation", &b5.signal_fluctuation, "signal_fluctuation/O");
   tree5->Branch("lpm_effect", &b5.lpm_effect, "lpm_effect/O");
   tree5->Branch("tau_regeneration", &b5.tau_regeneration, "tau_regeneration/O");
   tree5->Branch("planewaves", &b5.planewaves, "planewaves/O");
   tree5->Branch("depthdependetnfirn", &b5.depthdependentfirn, "depthdependentfirn/O");
   tree5->Branch("firn_presence", &b5.firn_presence, "firn_presence/O");
   tree5->Branch("spectrum", &b5.spectrum, "spectrum/O");
   tree5->Branch("widespectrum", &b5.widespectrum, "widespectrum/O");
   tree5->Branch("gzk_flux", &b5.gzk_flux, "gzk_flux/O");
   tree5->Branch("seckel_effect", &b5.seckel_effect, "seckel_effect/O");
   tree5->Branch("atten_freq", &b5.atten_freq, "atten_freq/I");
   tree5->Branch("atten_exp", &b5.atten_exp, "atten_exp/I");
   tree5->Branch("const_attenlength", &b5.const_attenlength, "const_attenlength/I");
   tree5->Branch("shadowing_dir_ref", &b5.shadowing_dir_ref, "shadowing_dir_ref/I");

//b5.nnu = NNU; //actually modified downstairs due to increase of stations
   b5.atgap = ATGap; //make sure to change to float if required
   b5.nrows = NROWS;
   b5.ncols = NCOLS;

   b5.signal_fluctuation = SIGNAL_FLUCT;
   b5.lpm_effect = LPM;
   b5.tau_regeneration = TAUREGENERATION;
   b5.planewaves = PLANEWAVE;
   b5.depthdependentfirn = DEPTH_DEPENDENT_N;
   b5.firn_presence = FIRN;
   b5.spectrum = SPECTRUM;
   b5.widespectrum = WIDESPECTRUM;
   b5.gzk_flux = GZK;
   b5.seckel_effect = SECKEL;
   b5.atten_freq = ATTEN_FREQ;
   b5.atten_exp = ATTEN_EXP;
   b5.const_attenlength = CONST_ATTENLENGTH;
   b5.shadowing_dir_ref = SHADOWING;

//tree5->Fill(); //re-used after NNU downstairs
//-----------------------------------------------------------------------------------------

//COMMENTED OUT ON 5/11/11 for smaller output files

//------------------------------------------------------------------+
   TTree* tree7 = new TTree("PAMe", "a tree for Energy tracking variables"); //
//------------------------------------------------------------------+
   tree7->Branch("N_Ant_perST", &b1.N_Ant_perST, "N_Ant_perST/I");

   tree7->Branch("Evmmhz1m_max", &b7.Evmmhz1m_max, "Evmmhz1m_max/F");
   tree7->Branch("Edeltheta_em_max", &b7.Edeltheta_em_max, "Edeltheta_em_max/F");
   tree7->Branch("Edeltheta_had_max", &b7.Edeltheta_had_max, "Edeltheta_had_max/F");
   tree7->Branch("En1", &b7.En1, "En1/F");
   tree7->Branch("Eshowerlength", &b7.Eshowerlength, "Eshowerlength/F");
   tree7->Branch("Evmmhz", &b7.Evmmhz, "Evmmhz[N_Ant_perST][95]/F");
   tree7->Branch("Evmmhz_mirror", &b7.Evmmhz_mirror, "Evmmhz_mirror[N_Ant_perST][95]/F");
   tree7->Branch("Esum_vmmhz", &b7.Esum_vmmhz, "Esum_vmmhz[N_Ant_perST]/F");
   tree7->Branch("Esum_vmmhz_mirror", &b7.Esum_vmmhz_mirror, "Esum_vmmhz_mirror[N_Ant_perST]/F");
   tree7->Branch("Evmmhz1m_unattened", &b7.Evmmhz1m_unattened, "Evmmhz1m_unattened[95]/F");
   tree7->Branch("Evmmhz1m_unattened_mirror", &b7.Evmmhz1m_unattened_mirror, "Evmmhz1m_unattened_mirror[95]/F");
   tree7->Branch("Esumvmmhz1m_unattened", &b7.Esumvmmhz1m_unattened, "Esumvmmhz1m_unattened/F");
   tree7->Branch("Esumvmmhz1m_unattened_mirror", &b7.Esumvmmhz1m_unattened_mirror, "Esumvmmhz1m_unattened_mirror/F");
   tree7->Branch("Edeltheta_em", &b7.Edeltheta_em, "Edeltheta_em[95]/F");
   tree7->Branch("Edeltheta_had", &b7.Edeltheta_had, "Edeltheta_had[95]/F");
   tree7->Branch("Eemfrac", &b7.Eemfrac, "Eemfrac/F");
   tree7->Branch("Ehadfrac", &b7.Ehadfrac, "Ehadfrac/F");
   tree7->Branch("Evmmhz_taper", &b7.Evmmhz_taper, "Evmmhz_taper[N_Ant_perST][95]/F");
   tree7->Branch("Evmmhz_taper_mirror", &b7.Evmmhz_taper_mirror, "Evmmhz_taper_mirrror[N_Ant_perST][95]/F");
   tree7->Branch("Esum_vmmhz_taper", &b7.Esum_vmmhz_taper, "Esum_vmmhz_taper[N_Ant_perST]/F");
   tree7->Branch("Esum_vmmhz_taper_mirror", &b7.Esum_vmmhz_taper_mirror, "Esum_vmmhz_taper_mirror[N_Ant_perST]/F");
   tree7->Branch("Evmmhz_em", &b7.Evmmhz_em, "Evmmhz_em[N_Ant_perST][95]/F");
   tree7->Branch("Evmmhz_had", &b7.Evmmhz_had, "Evmmhz_had[N_Ant_perST][95]/F");
   tree7->Branch("Evmmhz_em_mirror", &b7.Evmmhz_em_mirror, "Evmmhz_em_mirror[N_Ant_perST][95]/F");
   tree7->Branch("Evmmhz_had_mirror", &b7.Evmmhz_had_mirror, "Evmmhz_had_mirror[N_Ant_perST][95]/F");


//SHIFTED TO b7 ON 07/17/11 from b1
   tree7->Branch("term_LPA_my", &b7.term_LPA_my, "term_LPA_my[N_Ant_perST][95]/F");
   tree7->Branch("heff_my", &b7.heff_my, "heff_my[95]/F");
   tree7->Branch("vmmhz_my", &b7.vmmhz_my, "vmmhz_my[N_Ant_perST][95]/F");
   tree7->Branch("freq_my", &b7.freq_my, "freq_my[95]/F");
   tree7->Branch("vmmhz_mirror_my", &b7.vmmhz_mirror_my, "vmmhz_mirror_my[N_Ant_perST][95]/F");
   tree7->Branch("term_LPA_mirror_my", &b7.term_LPA_mirror_my, "term_LPA_mirror_my[N_Ant_perST][95]/F");
   tree7->Branch("Efield_atposnu_my", &b7.Efield_atposnu_my, "Efield_atposnu_my[N_Ant_perST][95]/F");
   tree7->Branch("Efield_atposnu_mirror_my", &b7.Efield_atposnu_mirror_my, "Efield_atposnu_mirror_my[N_Ant_perST][95]/F");


//--------------------------------------------------------------------------


//cout<<"energy is 10^"<<EXPONENT<<" eV"<<endl;
   time_t current_time;
   time(&current_time);
//if seed is set to zero, then the random seed will be set to the time.  Otherwise, it will be the seed from the input file
      if (seed == 0){
   int seed = time(NULL);
   }
//seed=343245;
//seed = (int) atoi(argv[4]); //KD: to input these two parameters from command line
   Rand3.SetSeed(seed);

   //Rand3.SetSeed(444456);  //Rand3.SetSeed(343245);
   //Rand3.SetSeed(current_time); // CJR - never set seed to time!

//KD:   output<<"pnu=1e"<<EXPONENT<<endl;//KD: 7/9/10, I was using this before, now commented out on 7/23
   pnu = pow(10, EXPONENT);
   double L_TauDecay = GetDecayLength(EXPONENT);
   //cout<<pnu<<endl;
   //get shower strength according to David Seckel's scaling
   double em_shower_length = 1.;
   double had_shower_length = 1.;

   if (SECKEL == 1) { //if we use David Seckel's scaling formula for the shower
      //read in the parameters for VmMHz1m
      //to get shower length

      double logx[25];
      double logl[25];
      logx[0] = 0;
      logl[0] = 0;

      ifstream lem("emhad2.txt");
      for (int k = 1; k < 25; k++) {
         lem >> logx[k] >> logl[k];
         if (log10(pnu / 1E15) < logx[k]) {
            em_shower_length = pow(10, logl[k - 1] + (logl[k] - logl[k - 1]) * (log10(pnu / 1E15) - logx[k - 1]) / (logx[k] - logx[k - 1]));
            //   cout<<"em_shower_length"<<em_shower_length<<endl;
            break;
         }
      }
      double x_energy[31];
      double tem_l[31];
      x_energy[0] = 0.;
      tem_l[0] = 0.;
      ifstream lhad("had_em_ratio.txt");
      for (int i = 1; i < 31; i++) {
         lhad >> x_energy[i] >> tem_l[i];
         if (log10(pnu / 1E15) < x_energy[i]) {
            had_shower_length = tem_l[i - 1] + (tem_l[i] - tem_l[i - 1]) * (log10(pnu / 1E15) - x_energy[i - 1]) / (x_energy[i] - x_energy[i - 1]);
            break;
         }
      }
      em_shower_length /= 1.01732;
      //  cout<<"em_shower_length"<<em_shower_length<<endl;
      had_shower_length /= 1.01732; //normalization
   }//END of if(SECKEL)

   //properties of antennas
   int iRow;
   int iCol;//the postion of an antenna
   int posnu_iRow = 0;
   int posnu_iCol = 0; //interaction bin
   double ATCoordinate[3];
   double MirrorATCoordinate[3];

   double ST_Posi_x[NROWS][NCOLS];
   double ST_Posi_y[NROWS][NCOLS];
   double ST_Posi_z[NROWS][NCOLS];

   string nuflavor; //neutrino flavor
   string current; //CC or NC?
   double emfrac;
   double hadfrac;
   double elpm;//LPM energy
   double elast_y;//inelasticity

   //Amy's Gaussian distribution of signal
   double vmmhz1m_max; // maximum V/m/MHz at 1m from Jaime (highest frequency)
   double vmmhz_max;//maximum V/m/MHz,  E field
   double vmmhz_max_freq[NFREQ];//maximum E field as a function of frequency
   double deltheta_em[NFREQ], deltheta_had[NFREQ];
   double vmmhz[NFREQ]; //  V/m/MHz at balloon (after all steps)
   double vmmhz1m[NFREQ];//for Seckel's scaling
   double vmmhz1m_unattened[NFREQ];//the E field at 1m of posnu before attenuation
   double vmmhz1m_unattened_max[NFREQ];
   double deltheta_em_max, deltheta_had_max;    // maximum value of above array angular distribution
   double eshower_em, eshower_had, showerlength_em;  //KD
   double vmmhz_em[NFREQ];
   double vmmhz_had[NFREQ]; //in Taper function
   //signals
//  double hitangle_e_ST0,hitangle_h_ST0;//ST_TYPE=0
   double hitangle_e_LPA, hitangle_h_LPA; //ST_TYPE=4



   for (int i = 0; i < NROWS; i++) {
      for (int j = 0; j < NCOLS; j++) {

         if (HEXAGONAL) {
            if (abs(i - j) == 0 || abs(i - j) == 2 || abs(i - j) == 4) {
               continue;
            }
         }

         GetATLocation(i, j, ATCoordinate);
         outantposall << ATCoordinate[0] << " " << ATCoordinate[1] << " " << ATCoordinate[2] << endl;
//KD: output all created antennas

// add PAMSt variable here
         ST_Posi_x[i][j] = ATCoordinate[0];
         ST_Posi_y[i][j] = ATCoordinate[1];
         ST_Posi_z[i][j] = ATCoordinate[2];
         b3.ST_Posi_x[i][j] = ST_Posi_x[i][j];
         b3.ST_Posi_y[i][j] = ST_Posi_y[i][j];
         b3.ST_Posi_z[i][j] = ST_Posi_z[i][j];
         //tree8->Fill();
      }
   }

//tree3->Fill();

   //for (ST_TYPE==0)
//  double term_e_ST0,term_h_ST0;
//  double e_component_ST0, h_component_ST0;

   // term for ST_TYPE 4, LPA
   double e_component_LPA, h_component_LPA;
   double term_LPA;
   double term_LPA_e, term_LPA_h;


   /*  //KD: commented out 05/17/11
   for (int i=0;i<4;i++) {
     bwslice_min[i]=bwslice_center[i]-bwslice_width[i]/2; // get low edge of bandwidth slices
     bwslice_max[i]=bwslice_center[i]+bwslice_width[i]/2;  // get upper edge of bandwidth slices
   }

   for(int i=0;i<4;i++)
     {
       bwslice_vnoise[i]=GetNoise(bwslice_max[i]-bwslice_min[i]);
     }
   */


   printf("NFREQ=%d, FREQ_LOW=%g, FREQ_HIGH=%g, BW=%g, FREQ_BIN=%g\n",
          NFREQ, FREQ_LOW, FREQ_HIGH, BW, FREQ_BIN);

   for (int i = 0; i < NFREQ; i++) {
      freq[i] = FREQ_LOW + (FREQ_HIGH - FREQ_LOW) * (double)i / (double)NFREQ; //MHz, freq. of each bin.
      b1.freq_my[i] = freq[i];
      b7.freq_my[i] = freq[i];
   }



   double Veff = 0.; //km3sr
   double E2dNdE = 0; //flux_limit
   double weight = 0.;
   // double dtryingdirection=0;
   double count_events = 0;
   //interaction point

   //ray tracing

   //double changle_deg=changle*RAD2DEG;

   double viewangle_mirror = 0.;
   // int ibin[2];//the bin of the interaction point. ibin[0] is the bin on x axis, ibin[1] is the bin on y axis.
   double viewangle = 0.; //the angle from the antenna to the axis(nnu) assuming the interaction point as the angle point.

   //neutrino path
   double nnu[3];//unit direction vector of neutrino
   // double nEfield[3];//unit direction vector of polarization
   //double lat_in=0;//latitude of entrance
   //double lon_in=0;//longitude of entrance
   //double chord=0;//from the earth entrance to the interaction point



   double posnu[3];
   double n_pol[3];
   double n_Bfield[3];
   double n_pol_out[3]; //for GetFresnel
   double posnu2AT[3];
   double posnu2MirrorAT[3];
   double nposnu2AT[3], nposnu2MirrorAT[3];
   double entry[3];
   double bounce[3]; //KD 06/29/11, use only in REF?
   double chord_inICE;//the path of a neutrino in the ice
   double chord_inCRUST;
   double chord2_inICE;
   double sigma;//cross section
   double int_len_kgm2;
   double attenlength_down;
   double attenlength_down_freq;
   double attenlength_up;
   double attenlength_up_freq;
   double integ_weight = 0;
   double abs_time;//absolute time the signal uses to reach the antenna
   double abs_time_mirror;
   double hy1 = 0, hy2 = 0; //the real path in each layer for direct events, h1--path in ice, h2--path in firn
   double hy1_mirror = 0, hy2_mirror = 0; //the real path in each layer for reflected events
   //for dipole
   double wavelength[NFREQ];
   double heff_dipole[NFREQ];//effective height of the dipole antenna
   double FACTOR_ofLAMDA = C / NICE / 1.e6;
   double FACTOR_ofHEFF = sqrt(2 * Zr * NICE * gain_dipole / Z0 / 4 / PI);
   if (DIPOLE) {
      for (int i = 0; i < NFREQ; i++) {
         wavelength[i] = FACTOR_ofLAMDA / freq[i];
         heff_dipole[i] = wavelength[i] * FACTOR_ofHEFF;
      }
   }

   double theta2 = 0, theta_nposnu2AT = 0, theta_nposnu2MirrorAT = 0, theta2_mirror = 0; //initial incident angle and the final refraction angle
   double phi_nposnu2AT = 0;
   double phi_nposnu2MirrorAT = 0; //phi angles for the vector from the interaction point to the AT
   double nsignal_mirror_atAT[3], nsignal_atAT[3]; //unit vector of the refracted signal at the AT position

   double theta_nposnu2ST = 0, theta_nposnu2MirrorST = 0; //KD introduced for ST plane wave calcs
   double ST_abs_time = 0, ST_abs_time_mirror = 0; //KD introduced for dir and ref abs times to ST
   double d_posnu2ST = 0, d_posnu2MirrorST = 0; //KD added for PLANEWAVE calcs
   double posnu2ST[3], posnu2MirrorST[3]; //KD see alse ~lines 1173
   double phi_nposnu2ST = 0, phi_nposnu2MirrorST = 0;
   double nsignal_atST[3], nsignal_atMirrorST[3]; //KD

   GetMaxDis(pnu, Max_distance);
   Max_distance = 4000.; //KD Max_distance=10000.;
   EDGE = Max_distance / 2.;
   double VOLUME = ICETHICK * pow((NCOLS - 1) * 300 + 2 * EDGE, 2.); //used as a scaling with 300m separation as benchmark, see NNU below

   double volume;
   if (HEXAGONAL)
// volume=ICETHICK*((NCOLS-1)*500+2*EDGE)*((NROWS-1)*866+2*EDGE); //total volume m^3
      volume = ICETHICK * ((2) * (1000.*HRAfactor) + 2 * EDGE) * ((2) * (1000.*HRAfactor) + 2 * EDGE);

   else
      volume = ICETHICK * pow((NCOLS - 1) * ATGap + 2 * EDGE, 2.); //total volume m^3

//  NNU=(int)(NNU*volume/VOLUME);
// 05/09/2011 KD: leaving NNU as absolute given number for the moment for multi-station


   //set ice thickness as 624m and set the distance between antennas as 300m
   //set 100*100 antennas  array

   //these two lines included here from above tree5
   b5.nnu = NNU;
   tree5->Fill();


   double TIMEPRECISION;
   TIMEPRECISION = 5.E-9; //setting uniform 5ns precision

   double VNOISE;
   VNOISE = GetNoise(BW);
   double NV = VNOISE * NSIGMA;

   GetBeamWidths(flare, gain, freq);


   printf("VNOISE=%g, NSIGMA=%g, NV=%g\n",
          VNOISE, NSIGMA, NV);

   printf("GAINFILENAME=[%s]\n", GAINFILENAME.c_str());

   ReadGains();

   // output<<"NNU="<<NNU<<endl;
   // output<<"ATGap="<<ATGap<<endl;
   // output<<"NFIRN="<<NFIRN<<endl;

   double n2;//index at the surface
   n2 = FIRN ? NFIRN : NICE; //KD: conditional operator means that if FIRN is present take NFIRN, otherwise take NICE

   int direct_events = 0, reflect_events = 0, combo_events = 0; //(KD:rootclean_A),both_events=0;
//KD:rootclean_A    double dir_w=0.,ref_w=0.,combo_w=0.,both_w=0.;
   double fan_dir_w = 0, fan_ref_w = 0, fan_combo_w = 0; //KD

//KD: added for flavor tracking
   int  events_nue = 0, events_numu = 0, events_nutau = 0;
   double  nue_w = 0, numu_w = 0, nutau_w = 0;
   double Veff_nue = 0., Veff_numu = 0., Veff_nutau = 0.;
   int nue_counts = 0, numu_counts = 0, nutau_counts = 0;

   //KD: I understand the following: check if it's within a reasonable meaningful viewangle of station
   //for the trigger level one(L1)---being on the cone that signals are possibly get there
   vector<int> iRow_oncone;
   vector<int> iCol_oncone;
   vector<int> iRow_oncone_mirror;
   vector<int> iCol_oncone_mirror;
   //the station get double hits, both direct signals and reflected signals can trigger this station
   vector<int> STtriggered_byboth_irow;
   vector<int> STtriggered_byboth_icol;

//KD: cout<<"SPECTRUM "<< SPECTRUM << " WIDESPECTRUM "<< WIDESPECTRUM << " GZK "<< GZK << endl;


//  int count_step=0; //KD:goes with definition just below


   for (int inu = 0; inu < NNU; inu++) {
      
      if ( (inu%10000)==0 ) {
         fprintf(stderr, "Processing %d / %d...            \r", inu, NNU);
      }
      

      if (SPECTRUM) {
         if (WIDESPECTRUM) {
            if (GZK)
               GetGZK();
            else
               GetE2();

            pnu = PickEnergy(energy, EdNdEdAdt); // KD : seems to work now, needed to modify array sizes in declaration.hh
         }


         else {
            if (GZK) {
               if (FANFLUX)
                  GetGZK_fan();
               else
                  GetGZK7();
            }

            else
               GetE2_fan();

            pnu = PickEnergy_fan(energy, EdNdEdAdt);

         }

         EXPONENT = log10(pnu);

      } //END of if(SPECTRUM)



      Zero(posnu, 3);
      Zero(nnu, 3);
      Zero(entry, 3);

//posnu[0]=(100. + 300.*(-2*Rand3.Rndm()+1)) ; posnu[1]=0. ; posnu[2]=300.  ;
      GetInteractionPoint(posnu);

//posnu[0]=342.098 ; posnu[1]= -142.388; posnu[2]= 507.861;//KD fringe neutrino example 1
//posnu[0]=240.364 ; posnu[1]= 870.649; posnu[2]= 127.351;  //KD neutrino example 2

      //GetPosnuBin(posnu[0],posnu[1],ibin);
      // test<<ibin[0]<<"  1  "<<ibin[1]<<endl;


//KD: added to take firn into account, imposed condition on 01/19/11
      if (FIRN) {
         double mytemp = Rand3.Rndm();
         if (posnu[2] > (ICETHICK - FIRNDEPTH) && mytemp > 0.59782)
            //how does this get affected for DEPTH_DEPEND? we are now using some average
            //KD 02/17/2011: it's the weight factor correction because fewer nucleons in the firn. Now, for uniform firn (ie no DEPTH_DEPEND), that cut should be reflecting the uniform density of firn at 0.3gcm^-3, so strictly, this can be further refined to be depth dependent
            continue;
      }


//+ outposnuall<< posnu[0]<<" "<<posnu[1]<<" "<<posnu[2]<<endl;   //KD: output all generated neutrinos?

      GetAttenlength(posnu, attenlength_up, attenlength_down);
      // cout<<attenlength_up<<" "<<attenlength_down<<endl;

      n1 = GetN(posnu); //get the refraction index at the interaction point
      elpm = GetLPM(n1); //KD: probably shouldn't be n1 but apparently it gets overriden in GetSpread below anyway?
      // cout<<"posnu[2]="<<posnu[2]<<"  n1="<<n1<<endl;
      changle = acos(1 / n1);
      f0 = 2.53E-7 / freq0 / sin(changle); //f0 in VmMHz1m function

      //Get the direction of a neutrino by random
      double theta_nu = acos(-2 * Rand3.Rndm() + 1); //original over all sky 12/9/2010
      //double theta_nu=acos(-1*Rand3.Rndm()); //original over half sky


      //double theta_nu=acos((-2*Rand3.Rndm()+1)/20.);
      // h1angle->Fill(theta_nu);
      // h1cos->Fill(cos(theta_nu));
      double phi_nu = Rand3.Rndm() * 2 * PI;
//double phi_nu=1.25 + 0.05*(-2*Rand3.Rndm()+1);
//double phi_nu = 0;


      nnu[0] = sin(theta_nu) * cos(phi_nu);
      nnu[1] = sin(theta_nu) * sin(phi_nu);
      nnu[2] = cos(theta_nu);
//nnu[0]=-0.68131 ; nnu[1]= 0.710366;  nnu[2]=-0.176627; // KD fringe example neutrino 1
//nnu[0]=-0.784983   ; nnu[1]= -0.619009;    nnu[2]=-0.0250813; //KD neutrino example 2

//SphAngles(nnu,sphangles);
//cout<<"KD2a: nnu=["<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<"]"<<"; nnu sphangles: "<<sphangles[0]<<"rad, "<<sphangles[1]<<"rad = "<<sphangles[0]*RAD2DEG<<"deg, "<<sphangles[1]*RAD2DEG<<"deg"<<endl;

      // if(theta_nu>=PI/2)//study the downward neutrinos
      //neglecting upgoing neutrinos
      // continue;

      GetEntryPoint(theta_nu, nnu, posnu, entry);
      // outVector(entry,3);

      chord_inICE = Distance(posnu, entry); //path in the ice
      b1.entrypt[0] = entry[0];
      b1.entrypt[1] = entry[1];
      b1.entrypt[2] = entry[2];

      // path length in ice available for tau decay
      double L_Ice = posnu[2]/cos(theta_nu);
//+++++++++++++++++++++++++++++++++++++++++++++++++

      GetNuFlavor(nuflavor); //nuflavor="nutau";//nuflavor="nue";
      if (nuflavor == "nue") {
         nue_counts++; //KD
      } else if (nuflavor == "numu") {
         numu_counts++; //KD
      } else {
         nutau_counts++; //KD
      }



      GetCurrent(current); //current="cc";

      //elast_y=0.1;//
      elast_y = Gety(); //elast_y=0.999; //KD setting to a fixed value to see

      GetEmHadFrac(nuflavor, current, elast_y, theta_nu, L_TauDecay, L_Ice, emfrac, hadfrac);
//cout<<"KD1a: "<<"emfrac="<<emfrac<<"  hadfrac="<<hadfrac<<" elast_y="<<elast_y<<endl; //KD, was previously commented, 7/9/10

//+++++++++++++++++++++++++++++++++++++++++++++++++

      if (theta_nu <= PI / 2) { //upward neutrinos
         GetChord(theta_nu, chord_inCRUST, chord2_inICE);
      }

      Sigma(pnu, sigma, int_len_kgm2);

      //output<<"sigma="<<sigma<<endl;
      // cout<<sigma*DensityICE/AMU<<endl;
      //for testing the normalization of nnu
      //   double normal=Mag(nnu);
      // cout<<normal<<endl;
      // cout<<CLOCKS_PER_SEC<<endl;

      GetPosnuBin(posnu, posnu_iRow, posnu_iCol);
      int iRow_min, iRow_max, iCol_min, iCol_max; //the possible block where antennas might be triggered
      iRow_min = posnu_iRow - (int)(Max_distance / ATGap);
      if (iRow_min < 0)
         iRow_min = 0;
      iRow_max = posnu_iRow + (int)(Max_distance / ATGap);
      if (iRow_max > NROWS)
         iRow_max = NROWS;

      iCol_min = posnu_iCol - (int)(Max_distance / ATGap);
      iCol_max = posnu_iCol + (int)(Max_distance / ATGap);
      if (iCol_min < 0)
         iCol_min = 0;
      if (iCol_max > NCOLS)
         iCol_max = NCOLS;

      iRow_oncone.clear();
      iCol_oncone.clear();
      iRow_oncone_mirror.clear();
      iCol_oncone_mirror.clear();
      STtriggered_byboth_irow.clear();
      STtriggered_byboth_icol.clear();


      for (iRow = iRow_min; iRow < iRow_max; iRow++) {
         for (iCol = iCol_min; iCol < iCol_max; iCol++) {

            if (HEXAGONAL) {
               //KD making sure we skip these alternate stations in hexagonal setup
               //it doesn't quite matter that we use ATGap of 500 for x and y above.
               if (abs(iRow - iCol) == 0 || abs(iRow - iCol) == 2 || abs(iRow - iCol) == 4)
                  continue;
            }

            GetBothLocation(iRow, iCol, ATCoordinate, MirrorATCoordinate);

            if (Distance(ATCoordinate, posnu) > Max_distance)
               continue;

            VectorMinus(ATCoordinate, posnu, posnu2AT);


//KD: INTRODUCING PLANE WAVE CALCULATIONS 01/14/10
            if (PLANEWAVE) { //for direct
//copied from AT calculations
               d_posnu2ST = Mag(posnu2AT); //non-normalized vector to STcenter
               nVector(posnu2AT, posnu2ST); //normalized vector to STcenter
               SphAngles(posnu2ST, sphangles); //converting to polar coords
               phi_nposnu2ST = sphangles[0] ;//rad
               theta_nposnu2ST = sphangles[1];//rad
               b1.phi_nposnu2ST = phi_nposnu2ST * RAD2DEG;
               b1.theta_nposnu2ST = theta_nposnu2ST * RAD2DEG; //WARNING: IT GETS RE-EVALUATED IN CASE OF FIRN


// Those output can be activated if a direction to center of station needed
// Also create a mirror one
//b1.nposnu2ST_i=  posnu2ST[0];
//b1.nposnu2ST_j=  posnu2ST[1];
//b1.nposnu2ST_k=  posnu2ST[2];

               theta2 = theta_nposnu2ST;

               double theta1 = theta2;
               hy1 = d_posnu2ST;
               hy2 = 0.; //if no firn

//KD: 01/18/11 added these 3 lines for cases of NO FIRN but for PLANEWAVE below
               nsignal_atST[0] = sin(theta2) * cos(phi_nposnu2ST);
               nsignal_atST[1] = sin(theta2) * sin(phi_nposnu2ST);
               nsignal_atST[2] = cos(theta2);


               if (FIRN) {

                  if (posnu[2] >= (ICETHICK - FIRNDEPTH)) { //if direct events interact in the firn
                     hy2 = hy1; //distance from the interaction point to the STATION
                     hy1 = 0; //nothing in bulk ice, KD
                  }

                  else {
                     double x1 = 1.e-100;
                     double x3 = 0.;
                     double h1 = ICETHICK - FIRNDEPTH - posnu[2];
                     double h2 = FIRNDEPTH;
                     double nconst = n2 * n2 / n1 / n1; //here n2=NFIRN, n1=NICE
                     double x2 = nconst; // x2 is not angle. x2=sin(theta)^2
                     double deltax = sqrt(Square(posnu[0] - ATCoordinate[0]) + Square(posnu[1] - ATCoordinate[1])); //KD only xyplane distance

                     if (deltax == 0) {
                        theta1 = 0.;
                        theta2 = 0.;
                        hy1 = h1;
                        hy2 = h2;
                     }

                     else {
                        do {
                           x3 = (x1 + x2) / 2.;
                           if (fx(x1, h1, h2, nconst, deltax)*fx(x3, h1, h2, nconst, deltax) < 0)
                              x2 = x3;
                           else x1 = x3;
                        }

                        while (fabs(fx(x3, h1, h2, nconst, deltax)) > 0.00001);
                        theta1 = asin(sqrt(x3));
                        theta2 = asin(n1 * sin(theta1) / n2);
                        hy1 = h1 / cos(theta1);
                        hy2 = h2 / cos(theta2);
                     }

                     theta_nposnu2ST = theta1; //theta of the signal close to the interaction point
                     d_posnu2ST = hy1 + hy2;
                  }

                  nsignal_atST[0] = sin(theta2) * cos(phi_nposnu2ST);
                  nsignal_atST[1] = sin(theta2) * sin(phi_nposnu2ST);
                  nsignal_atST[2] = cos(theta2);

               }//end of if FIRN

               ST_abs_time = (hy1 * NICE + hy2 * NFIRN) / C; //entry in first vector is for direct
            } //end of if PLANEWAVE (for direct)

            b1.theta_signal_atAT = theta2 * RAD2DEG;


            //viewangle=Angle(nsignal,nnu);
            viewangle = Angle(posnu2AT, nnu);
//if(fabs(viewangle-changle)>=0.0)
            if (fabs(viewangle - changle) < CONEWIDTH) {
               iRow_oncone.push_back(iRow);//read out the information of first_step trigged antennas
               iCol_oncone.push_back(iCol);//here the count_Atoncone begins from 1 to the biggest number of the trigged number
            }



//QUICK CHECKS FOR MIRROR NOW
            //GetMirrorATLocation(iRow,iCol,MirrorATCoordinate);
            if (Distance(MirrorATCoordinate, posnu) > Max_distance)
               continue;

            // h_distance=sqrt(Square(posnu[0]-MirrorATCoordinate[0])+Square(posnu[1]-MirrorATCoordinate[1]));
            //need to get the incident angle theta1_mirror
            // theta2_mirror=asin(n1*sin(theta1_mirror)/n2);
            VectorMinus(MirrorATCoordinate, posnu, posnu2MirrorAT);

            if (PLANEWAVE) { //now for mirror
               d_posnu2MirrorST = Mag(posnu2MirrorAT);
               nVector(posnu2MirrorAT, posnu2MirrorST);
               SphAngles(posnu2MirrorST, sphangles);
               phi_nposnu2MirrorST = sphangles[0] ;//rad
               theta_nposnu2MirrorST = sphangles[1];//rad
               b1.phi_nposnu2MirrorST = phi_nposnu2MirrorST * RAD2DEG;
               b1.theta_nposnu2MirrorST = theta_nposnu2MirrorST * RAD2DEG; //WARNING: IT GETS RE-EVALUATED IN CASE OF FIRN

               theta2 = theta_nposnu2MirrorST;

               double theta1_mirror = theta2_mirror;
               hy1_mirror = d_posnu2MirrorST;
               hy2_mirror = 0.; //if no firn
//KD: 01/18/11 added these 3 lines for cases of NO FIRN but for PLANEWAVE below
               nsignal_atMirrorST[0] = sin(theta2_mirror) * cos(phi_nposnu2MirrorST);
               nsignal_atMirrorST[1] = sin(theta2_mirror) * sin(phi_nposnu2MirrorST);
               nsignal_atMirrorST[2] = cos(theta2_mirror);

               if (FIRN) {

                  if (posnu[2] < (ICETHICK - FIRNDEPTH)) { //interact in the ice
                     double x1_mirror = 1.e-100;
                     double x3_mirror = 0.;
                     double h1_mirror = (ICETHICK - FIRNDEPTH) + posnu[2];
                     double h2_mirror = FIRNDEPTH;
                     double nconst_mirror = n2 * n2 / n1 / n1; //here n2=NFIRN, n1=NICE
                     double x2_mirror = nconst_mirror;
                     double deltax_mirror = sqrt(Square(posnu[0] - MirrorATCoordinate[0]) + Square(posnu[1] - MirrorATCoordinate[1]));

                     if (deltax_mirror == 0) {
                        theta1_mirror = 0.;
                        theta2_mirror = 0.;
                        hy1_mirror = h1_mirror;
                        hy2_mirror = h2_mirror;
                     } else {
                        do {
                           x3_mirror = (x1_mirror + x2_mirror) / 2;
                           if (fx(x1_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)*fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror) < 0)
                              x2_mirror = x3_mirror;
                           else x1_mirror = x3_mirror;
                        } while (fabs(fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)) > 0.0001);
                        theta1_mirror = asin(sqrt(x3_mirror)); //incident angle in the ice,less than 90 degrees
                        theta2_mirror = asin(n1 * sin(theta1_mirror) / n2); //refraction angle in the firn, less than 90 degrees
                        hy1_mirror = h1_mirror / cos(theta1_mirror);
                        hy2_mirror = h2_mirror / cos(theta2_mirror);
                     }

                     theta_nposnu2MirrorST = PI - theta1_mirror; //interacting in the ice, the theta angle close to the interaction point
                     d_posnu2MirrorST = hy1_mirror + hy2_mirror;
                  }//end of if(<400m), interact in the ice

                  else { //reflected events happen in the firn
                     double x1_mirror = 1.e-100;
                     double x3_mirror = 0.;
                     double h1_mirror = 2 * (ICETHICK - FIRNDEPTH);
                     double h2_mirror = FIRNDEPTH + (posnu[2] - (ICETHICK - FIRNDEPTH));
                     //double nconst_mirror=n2*n2/n1/n1;//this is wrong because in this case n1=n2=NFIRN since the interaction point is in the firn
                     double nconst_mirror = NFIRN * NFIRN / NICE / NICE; //for FIRN case
                     double x2_mirror = nconst_mirror;
                     double deltax_mirror = sqrt(Square(posnu[0] - MirrorATCoordinate[0]) + Square(posnu[1] - MirrorATCoordinate[1]));
                     if (deltax_mirror == 0) {
                        theta1_mirror = 0.;
                        theta2_mirror = 0.;
                        hy1_mirror = h1_mirror;
                        hy2_mirror = h2_mirror;
                     } else {
                        do {
                           x3_mirror = (x1_mirror + x2_mirror) / 2;
                           if (fx(x1_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)*fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror) < 0)
                              x2_mirror = x3_mirror;
                           else x1_mirror = x3_mirror;
                        } while (fabs(fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)) > 0.00001);
                        theta1_mirror = asin(sqrt(x3_mirror)); //incident angle in the ice,less than 90 degrees
                        //theta2_mirror=asin(n1*sin(theta1_mirror)/n2);//refraction angle in the firn, less than 90 degrees
                        theta2_mirror = asin(NICE * sin(theta1_mirror) / NFIRN);
                        hy1_mirror = h1_mirror / cos(theta1_mirror);
                        hy2_mirror = h2_mirror / cos(theta2_mirror);
                     }
                     theta_nposnu2MirrorST = PI - theta2_mirror; //the theta angle at the mirror antennas, also the theta angle of the signal at the interaction point since interating in the firn
                     d_posnu2MirrorST = hy1_mirror + hy2_mirror;
                  }

//these might need to be double checked
                  nsignal_atMirrorST[0] = sin(theta2_mirror) * cos(phi_nposnu2MirrorST);
                  nsignal_atMirrorST[1] = sin(theta2_mirror) * sin(phi_nposnu2MirrorST);
                  nsignal_atMirrorST[2] = -cos(theta2_mirror);

               }// end of if FIRN

               ST_abs_time_mirror = (hy1_mirror * NICE + hy2_mirror * NFIRN) / C;
            }// end of if PLANEWAVE for reflected


            b1.theta_signal_atAT_mirror = theta2_mirror * RAD2DEG;


            // viewangle_mirror=Angle(nsignal_mirror,nnu);
            viewangle_mirror = Angle(posnu2MirrorAT, nnu);
            //if(fabs(viewangle_mirror-changle)>=0.0)             //this goes through all stations
            if (fabs(viewangle_mirror - changle) < CONEWIDTH) {
               iRow_oncone_mirror.push_back(iRow);
               iCol_oncone_mirror.push_back(iCol);
            }

         }//end of icol loop
      }//end of irow loop

      int sum_AToncone = iRow_oncone.size() + iRow_oncone_mirror.size();
      //  if((int)iRow_oncone_mirror.size()==0&&(int)iRow_oncone.size()==0)  {outVector(nnu,3); outVector(posnu,3);}
      // cout<<sum_AToncone<<endl;


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// KD What does this do?
//      if (ST_TYPE == 4) {
         if (sum_AToncone == 0) //KD: for ST_TYPE 4, looks like we just want more than one triggering?
//cout<<"KD : This is printed if sum_AToncone=0."<<endl;   //THIS COMMENT TERMINATES THE LOOP FOR SOME REASON
            continue;
//      }


      /* //KD: deactivated on 5/17/11. Still don't understand why should be >=3
      else
      {if(sum_AToncone<3) //KD: for others, why at least three?
      //cout<<"KD : sum_AToncone is between 0 and 3. Relevant for multistations?"<<endl;
        continue;
      }
      */
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      //begin to analyze the signals
      //first get the maximum signal and the spreading angles
      vmmhz1m_max = GetVmMHz1m(pnu, FREQ_HIGH, X0ICE, ECICE, n1, AEX_ICE);
      GetSpread(pnu, emfrac, hadfrac, FREQ_LOW, n1, X0ICE, deltheta_em_max, deltheta_had_max, eshower_em, eshower_had, showerlength_em); //maximum angular deviation from Cherenkov cone where signal is still detectable. //KD: deltheta are in radians

//-------------------------------
      b7.Evmmhz1m_max = vmmhz1m_max; //actually this is redundant, it's already in the file
      b7.Edeltheta_em_max = deltheta_em_max;
      b7.Edeltheta_had_max = deltheta_had_max;
      b7.En1 = n1;
      b7.Eshowerlength = showerlength_em;
//-------------------------------


      if (SECKEL == 0) {
         for (int i = 0; i < NFREQ; i++) {
            GetSpread(pnu, emfrac, hadfrac, freq[i], n1, X0ICE, deltheta_em[i], deltheta_had[i], eshower_em, eshower_had, showerlength_em);
//-------------------------------
            b7.Edeltheta_em[i]  = deltheta_em[i];
            b7.Edeltheta_had[i] = deltheta_had[i];
//-------------------------------
         }
      }
//-------------------------------
      b7.Eemfrac = emfrac;
      b7.Ehadfrac = hadfrac;
//-------------------------------


      int sum_triggeredST = 0;
      int sum_triggeredST_mirror = 0;
      int all_triggeredST = 0;
      int sum_ST = 0; //KD: relevant for Jiwoo's b1 , and only for ST_TYPE<2? on 01/13/11, KD modified it for PLANEWAVE

      int iAT_index_mymax = 0; //KD 01/28/11 to track antenna with max voltage
      int iAT_mirror_index_mymax = 0; //KD 01/28/11 to track antenna with max voltage
      int iAT_index_mymin = 0;
      int iAT_mirror_index_mymin = 0;
      int iAT_index_timemax = 0; //KD to track antenna with latest time
      int iAT_mirror_index_timemax = 0; //KD to track antenna with latest time
      int iAT_index_timemin = 0; //KD to track antenna with earliest time
      int iAT_mirror_index_timemin = 0; //KD to track antenna with earliest time

      //for plotting
      vector<double> viewangle_triggered;
      vector<double> viewangle_triggered_mirror;
      vector<double> attenfactor;
      vector<double> attenfactor_mirror;
      vector<double> dis;
      vector<double> dis_mirror;
      vector<double> vdis;
      vector<double> vdis_mirror;
      vector<double> cosz;
      vector<double> cosz_mirror;
      vector<double> coszp;
      vector<double> coszp_mirror;

      vector<double> theta_signal;
      vector<double> theta_signal_atAT;
      vector<double> theta_signal_mirror;
      vector<double> theta_signal_atAT_mirror;

      double new_posnu2AT[3];//after filtered by L1 trigger(whether on cone)
      double d_posnu2AT;
      double new_posnu2MirrorAT[3];
      double d_posnu2MirrorAT;



//KD: N_Ant_perST now overriden by input file. It will only affect ST_TYPE 4 as it stands.
//      int N_Ant_perST=1;//number of antennas each station
//      if (ST_TYPE == 3)
//         N_Ant_perST = 5;
//     if(ST_TYPE==4)
// N_Ant_perST=4;//to be changed accordingly //FW:8 log periodic antennas per station

//      if (ST_TYPE == 3 || ST_TYPE == 4) {

         //KD: added 01/13/2011
         //vector<double> ST_Posi_x;
         //vector<double> ST_Posi_y;
         //vector<double> ST_Posi_z;
         //vector<int> ST_iRow;
         //vector<int> ST_iCol;
         //ST_Posi_x.clear();
         //ST_Posi_y.clear();
         //ST_Posi_z.clear();
         //ST_iRow.clear();
         //ST_iCol.clear();


         //variables after event trigger(L3),including both direct and reflected events
         vector<int> N_TriggeredAT_perST_L3;//the number of triggered antennas in each triggered station
         vector<int> Trig_Type_L3;//triggered by direct or reflected signal;
         vector<int> iRow_L3;
         vector<int> iCol_L3;
         vector<double> Emax_atposnu_L3;
         //into each triggered station, analyze the antennas in the station
         vector<int> iAT_L3;//in the triggered station, which antennas are triggered, this is the label of the triggered antenna
         vector<double> path_inice_L3;
         vector<double> path_infirn_L3;
         vector<double> totaltime_L3;//the time for a signal getting to the antenna
         vector<double> Efield_atposnu_L3;

         vector<double> volte_allAT_L3;
         vector<double> volth_allAT_L3;
         vector<double> AT_posi_x_L3;
         vector<double> AT_posi_y_L3;
         vector<double> AT_posi_z_L3;
         vector<double> Dir_shower_x_L3;
         vector<double> Dir_shower_y_L3;
         vector<double> Dir_shower_z_L3;
         vector<double> theta_shower_atposnu_L3;

         N_TriggeredAT_perST_L3.clear();
         Trig_Type_L3.clear();
         iRow_L3.clear();
         iCol_L3.clear();

         iAT_L3.clear();
         path_inice_L3.clear();
         path_infirn_L3.clear();
         totaltime_L3.clear();
         Efield_atposnu_L3.clear();
         Emax_atposnu_L3.clear();
         volte_allAT_L3.clear();
         volth_allAT_L3.clear();
         AT_posi_x_L3.clear();
         AT_posi_y_L3.clear();
         AT_posi_z_L3.clear();
         Dir_shower_x_L3.clear();
         Dir_shower_y_L3.clear();
         Dir_shower_z_L3.clear();
         theta_shower_atposnu_L3.clear();



         //variables after station trigger(L2)
         vector<int> N_TriggeredAT_perST_L2;//the number of triggered antennas in each triggered station
         vector<int> Trig_Type_L2;//triggered by direct or reflected signal;
         vector<int> iRow_L2;
         vector<int> iCol_L2;

         //into each triggered station, analyze the antennas in the station
         vector<int> iAT_L2;//in the triggered station, which antennas are triggered, this is the label of the triggered antenna
         vector<double> path_inice_L2;
         vector<double> path_infirn_L2;
         vector<double> totaltime_L2;//the time for a signal getting to the antenna
         vector<double> Efield_atposnu_L2;
         vector<double> volte_allAT_L2;
         vector<double> volth_allAT_L2;
         vector<double> AT_posi_x_L2;
         vector<double> AT_posi_y_L2;
         vector<double> AT_posi_z_L2;
         vector<double> Dir_shower_x_L2;
         vector<double> Dir_shower_y_L2;
         vector<double> Dir_shower_z_L2;
         vector<double> theta_shower_atposnu_L2;
         vector<double> hitangle_e_LPA_L2;
         vector<double> hitangle_h_LPA_L2;

         N_TriggeredAT_perST_L2.clear();
         Trig_Type_L2.clear();
         iRow_L2.clear();
         iCol_L2.clear();


         iAT_L2.clear();
         path_inice_L2.clear();
         path_infirn_L2.clear();
         totaltime_L2.clear();
         Efield_atposnu_L2.clear();
         volte_allAT_L2.clear();
         volth_allAT_L2.clear();
         AT_posi_x_L2.clear();
         AT_posi_y_L2.clear();
         AT_posi_z_L2.clear();
         Dir_shower_x_L2.clear();
         Dir_shower_y_L2.clear();
         Dir_shower_z_L2.clear();
         theta_shower_atposnu_L2.clear();
         hitangle_e_LPA_L2.clear();
         hitangle_h_LPA_L2.clear();

         sum_triggeredST = 0;
         for (uint WhichStation = 0; WhichStation < iRow_oncone.size(); WhichStation++) {
            GetATLocation(iRow_oncone.at(WhichStation), iCol_oncone.at(WhichStation), ATCoordinate); //get the coordinates of the center of the station

//-----------------------------------------------------------------------------------------------------------
// ESTABLISHING MAIN TRIGGER + LOOSE TRIGGER CONDITIONS FOR MULTI-STATION
//-----------------------------------------------------------------------------------------------------------
//cout<<inu<<"():"<<iRow_oncone.at(WhichStation)<<","<<iCol_oncone.at(WhichStation)<<endl; // print out to establish varying station row/col
            if (HEXAGONAL) {
               if (iRow_oncone.at(WhichStation) == 1 && iCol_oncone.at(WhichStation) == 2) {
                  NV = VNOISE * NSIGMA;
               } else {
                  NV = VNOISE * NSIGMA; //Note: JCH changes this from 5.0 to NSIGMA July 4th 2012, because I don't understand how to fully control the threshold
               }
            }

//cout<<inu<<"(DIR):"<<iRow_oncone.at(WhichStation)<<","<<iCol_oncone.at(WhichStation)<<","<<NV;
//-----------------------------------------------------------------------------------------------------------



            /*        // KD COMMENTED OUT 5/9/11
                     //if st_type=3
                     double  ATCoordinates5[5][3]={{ATCoordinate[0],ATCoordinate[1],ATCoordinate[2]},{ATCoordinate[0]+ST_DISTANCE, ATCoordinate[1],ATCoordinate[2]},{ATCoordinate[0],ATCoordinate[1]-ST_DISTANCE, ATCoordinate[2]},
                               {ATCoordinate[0]-ST_DISTANCE, ATCoordinate[1],ATCoordinate[2]},{ATCoordinate[0],ATCoordinate[1]+ST_DISTANCE, ATCoordinate[2]}};
            */

            //Set Antenna Positions
            double ATCoordinates8[N_Ant_perST][3];//the detailed position of the center of each LPA in a station                        
            if (StationType == 0){ //All antennas pointing down, equally spaced around station center
                for (int i = 0; i < N_Ant_perST; i++) {
                   double phi = (2. / N_Ant_perST) * PI * i; //the phi angle of each LPA's center
                   ATCoordinates8[i][0] = ATCoordinate[0] + ST4_R * cos(phi);
                   ATCoordinates8[i][1] = ATCoordinate[1] + ST4_R * sin(phi);
                   ATCoordinates8[i][2] = ATCoordinate[2];
                }
            }
            else if (StationType == 1){ //Some custom antenna config
                cout<<"Haven't defined this yet"<<endl;
            }
            else {
                cout<<"Invalid Station type "<<StationType<<endl;
            }
            
            //Set Antenna Types
            int AntType[N_Ant_perST];
            //type 0 = 100MHz theoretical LPDA (original ShelfMC model)
            //type 1 = 100MHz Create LPDA, Anna's WhippleD model
            if (StationType == 0){
                for (int i = 0; i < N_Ant_perST; i++) {
                    AntType[i]=0;
                }
            }
            else if (StationType == 1){
                for (int i = 0; i < N_Ant_perST; i++) {
                    AntType[i]=1;
                }
            }
            else {
                cout<<"Invalid Station type "<<StationType<<endl;
            }

            //Set Antenna Orientation
            double Ant_n_boresight[N_Ant_perST][3];
            double Ant_n_eplane[N_Ant_perST][3];

            if (StationType == 0) {
                for (int i = 0; i < N_Ant_perST; i++) {
                    Ant_n_eplane[i][0] = cos((0.5 + i * (2. / N_Ant_perST))*PI);
                    Ant_n_eplane[i][1] = sin((0.5 + i * (2. / N_Ant_perST))*PI);
                    Ant_n_eplane[i][2] = 0.;
                    Ant_n_boresight[i][0] = 0.;
                    Ant_n_boresight[i][1] = 0.;
                    Ant_n_boresight[i][2] = -1.;
                }
            }
            else if (StationType == 1) {
                cout<<"haven't defined this yet"<<endl;
            }
            else {
                cout<<"Invalid Station Type"<<endl;
            }

            //define some variables before going into the antenna loop of one station

            //variables after antenna trigger(L1)

            //into each triggered station, analyze the antennas in the station
            vector<int> iAT_L1;//in the triggered station, which antennas are triggered, this is the label of the triggered antenna
            vector<double> path_inice_L1;
            vector<double> path_infirn_L1;
            vector<double> totaltime_L1;//the time for a signal getting to the antenna
            vector<double> Efield_atposnu_L1;
            vector<double> volte_allAT_L1;
            vector<double> volth_allAT_L1;
            vector<double> AT_posi_x_L1;
            vector<double> AT_posi_y_L1;
            vector<double> AT_posi_z_L1;
            vector<double> Dir_shower_x_L1;
            vector<double> Dir_shower_y_L1;
            vector<double> Dir_shower_z_L1;
            vector<double> theta_shower_atposnu_L1;
            vector<double> hitangle_e_LPA_L1;
            vector<double> hitangle_h_LPA_L1;

            iAT_L1.clear();
            path_inice_L1.clear();
            path_infirn_L1.clear();
            totaltime_L1.clear();
            Efield_atposnu_L1.clear();
            volte_allAT_L1.clear();
            volth_allAT_L1.clear();
            AT_posi_x_L1.clear();
            AT_posi_y_L1.clear();
            AT_posi_z_L1.clear();
            Dir_shower_x_L1.clear();
            Dir_shower_y_L1.clear();
            Dir_shower_z_L1.clear();
            theta_shower_atposnu_L1.clear();
            hitangle_e_LPA_L1.clear();
            hitangle_h_LPA_L1.clear();

//----KD added 7/8
            //for plotting, !!!!!!!!!!!!!!check whether the mirror part needs to be here too, or do I need to copy them in the mirror section
            viewangle_triggered.clear();
            viewangle_triggered_mirror.clear();
            attenfactor.clear();
            attenfactor_mirror.clear();
            dis.clear();
            dis_mirror.clear();
            vdis.clear();
            vdis_mirror.clear();
            cosz.clear();
            cosz_mirror.clear();
            coszp.clear();
            coszp_mirror.clear();

            theta_signal.clear();
            theta_signal_atAT.clear();
            theta_signal_mirror.clear();
            theta_signal_atAT_mirror.clear();
//---from ST<3--------


            //start the antennas loop of one station
            int sum_triggeredAT = 0;
            for (int WhichAntenna = 0; WhichAntenna < N_Ant_perST; WhichAntenna++)
// KD above was used for general trigger, but I'm manually modifying for 2 out of 3 case
//       for(int WhichAntenna=0; WhichAntenna<3; WhichAntenna++)
            {
               /* //COMMENTED OUT 5/9/11
               if(ST_TYPE==3)
                 VectorMinus(ATCoordinates5[WhichAntenna],posnu,new_posnu2AT);
                 */
//               if (ST_TYPE == 4)
               VectorMinus(ATCoordinates8[WhichAntenna], posnu, new_posnu2AT); //cout<<"I'm here 1"<<endl;
               nVector(new_posnu2AT, nposnu2AT);

               if (nposnu2AT[1] >= 0) {
                  if (nposnu2AT[0] > 0)
                     phi_nposnu2AT = atan(nposnu2AT[1] / nposnu2AT[0]);
                  else if (nposnu2AT[0] < 0)
                     phi_nposnu2AT = PI + atan(nposnu2AT[1] / nposnu2AT[0]);
                  else
                     phi_nposnu2AT = PI / 2.;
               } else {
                  if (nposnu2AT[0] > 0)
                     phi_nposnu2AT = 2 * PI + atan(nposnu2AT[1] / nposnu2AT[0]);
                  else if (nposnu2AT[0] < 0)
                     phi_nposnu2AT = PI + atan(nposnu2AT[1] / nposnu2AT[0]);
                  else
                     phi_nposnu2AT = 3 * PI / 2.;
               }

               theta_nposnu2AT = acos(nposnu2AT[2]);
               d_posnu2AT = Mag(new_posnu2AT);
               viewangle = Angle(nnu, nposnu2AT);
               GetPolarization(nnu, nposnu2AT, n_pol, n_Bfield);
               theta2 = theta_nposnu2AT;

               /*
               if (WhichAntenna==7){
               SphAngles(nposnu2AT,sphangles);
               cout<<
               "inu: "<<inu<<" nnu:"<< nnu[0] <<","<<nnu[1]<<","<<nnu[2]<<
               " nposnu2AT= ["<< nposnu2AT[0]<<","<<nposnu2AT[1]<<","<<nposnu2AT[2]<<"] " <<sphangles[0]*RAD2DEG<<"deg, "<<sphangles[1]*RAD2DEG<<"deg"<<
               " theta_nposnu2AT: "<<theta_nposnu2AT*RAD2DEG<<" viewangle "<<viewangle*RAD2DEG<<endl;
               //KD
               SphAngles(n_pol,sphangles);
               cout<<"inu: "<<inu<<" "<<"n_pol[]=["<<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<"]"<<
               " n_pol "<<sphangles[0]*RAD2DEG<<"deg, "<<sphangles[1]*RAD2DEG<<"deg"<<endl;
               }
               */
               double theta1 = theta2;
               hy1 = d_posnu2AT;
               hy2 = 0.; //if no firn

               //KD: 01/18/11 added these 3 lines for cases of NO FIRN but for PLANEWAVE below
               nsignal_atAT[0] = sin(theta2) * cos(phi_nposnu2AT);
               nsignal_atAT[1] = sin(theta2) * sin(phi_nposnu2AT);
               nsignal_atAT[2] = cos(theta2);

               if (FIRN) {

                  if (posnu[2] >= (ICETHICK - FIRNDEPTH)) { //if direct events interact in the firn
                     hy2 = hy1; //the distance from the interaction point to the antenna
                     hy1 = 0;
                  }

                  else { //if direct events interact in the ice
                     double x1 = 1.e-100;
                     double x3 = 0.;
                     double h1 = ICETHICK - FIRNDEPTH - posnu[2];
                     double h2 = FIRNDEPTH;
                     double nconst = n2 * n2 / n1 / n1; //here n2=NFIRN, n1=NICE
                     double x2 = nconst; // x2 is not angle. x2=sin(theta)^2
                     //double deltax=sqrt(Square(posnu[0]-ATCoordinate[0])+Square(posnu[1]-ATCoordinate[1]));
                     double deltax;
//                     if (ST_TYPE == 4)
                        deltax = sqrt(Square(posnu[0] - ATCoordinates8[WhichAntenna][0]) + Square(posnu[1] - ATCoordinates8[WhichAntenna][1]));
                     /* //COMMENTED OUT 5/9/11
                     else if(ST_TYPE==3)
                       deltax=sqrt(Square(posnu[0]-ATCoordinates5[WhichAntenna][0])+Square(posnu[1]-ATCoordinates5[WhichAntenna][1]));
                       */
//                     else {
//                        deltax = 0;
//                        cout << "Wrong ST_TYPE value" << endl;
//                     }
                     if (deltax == 0) {
                        theta1 = 0.;
                        theta2 = 0.;
                        hy1 = h1;
                        hy2 = h2;
                     } else {
                        do {
                           x3 = (x1 + x2) / 2.;

                           if (fx(x1, h1, h2, nconst, deltax)*fx(x3, h1, h2, nconst, deltax) < 0)
                              x2 = x3;
                           else x1 = x3;
                           // ncount++;
                           // if(ncount > 100) cout<<"ncount="<<ncount<<endl;
                        } while (fabs(fx(x3, h1, h2, nconst, deltax)) > 0.00001);

                        theta1 = asin(sqrt(x3));
                        theta2 = asin(n1 * sin(theta1) / n2);
                        hy1 = h1 / cos(theta1);
                        hy2 = h2 / cos(theta2);

                     }

                     theta_nposnu2AT = theta1; //theta of the signal close to the interaction point

                     d_posnu2AT = hy1 + hy2;

                     double original_nsignal[3];
                     original_nsignal[0] = sin(theta1) * cos(phi_nposnu2AT);
                     original_nsignal[1] = sin(theta1) * sin(phi_nposnu2AT);
                     original_nsignal[2] = cos(theta1);
                     viewangle = Angle(nnu, original_nsignal);

                     //KD introduced here 07/26/11 to get original modified polarization
                     GetPolarization(nnu, original_nsignal, n_pol, n_Bfield);
                     b1.original_nsignal[0] = original_nsignal[0];
                     b1.original_nsignal[1] = original_nsignal[1];
                     b1.original_nsignal[2] = original_nsignal[2];
                     b1.Polarization_D_original[0] = n_pol[0];
                     b1.Polarization_D_original[1] = n_pol[1];
                     b1.Polarization_D_original[2] = n_pol[2];
                     b1.Bfield_D_original[0] = n_Bfield[0];
                     b1.Bfield_D_original[1] = n_Bfield[1];
                     b1.Bfield_D_original[2] = n_Bfield[2];
                     SphAngles(n_pol, sphangles);
                     b1.polariz_phi_D_original = sphangles[0] * RAD2DEG;
                     b1.polariz_theta_D_original = sphangles[1] * RAD2DEG;

                     nsignal_atAT[0] = sin(theta2) * cos(phi_nposnu2AT);
                     nsignal_atAT[1] = sin(theta2) * sin(phi_nposnu2AT);
                     nsignal_atAT[2] = cos(theta2);

                     /*
                     cout<<" >after firn "<<inu<<" nnu:"<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P2:"
                        <<nsignal_atAT[0]<<","<<nsignal_atAT[1]<<","<<nsignal_atAT[2]<<" angles:"
                        <<theta2*RAD2DEG<<","<<phi_nposnu2AT*RAD2DEG<<" E1 tree:"
                        <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]
                        <<endl;
                     */


                     if (DEPTH_DEPENDENT_N) {
                        GetRotation(theta1, theta2, original_nsignal, n_pol, n_pol_out);

                        /*
                        cout<<" >after rot "<<inu<<" nnu:"<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P1:"
                           <<original_nsignal[0]<<","<<original_nsignal[1]<<","<<original_nsignal[2]<<" Snells angles:"
                           <<theta2*RAD2DEG<<"<->"<<theta1*RAD2DEG<<" E1:"
                           <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<" E2:"
                           <<n_pol_out[0]<<","<<n_pol_out[1]<<","<<n_pol_out[2]
                           <<endl;
                        */

                     } else {
//GetRotation(theta1, theta2, original_nsignal, n_pol, n_pol_out);
//use above if we want to use rotation for binary indices

//KD, now using what I think is correct way to get Fresnel'd polarization
// Is nsignal_atAT normalized?
                        //double n_pol_out[3];
                        GetFresnel(theta1, theta2, original_nsignal, nsignal_atAT, n_pol, n_pol_out);
                        //cout<<"n_pol B "<<n_pol[0]<<" "<<n_pol[1]<<" "<<n_pol[2]<<endl;
                        //cout<<"n_pol out B "<<n_pol_out[0]<<" "<<n_pol_out[1]<<" "<<n_pol_out[2]<<endl;
                        //cout<<"theta's B "<<theta1<<" "<<theta2<<endl;
                     }//end of DEPTH_DEPENDENT_N


// we feed the n_pol_out into the tree; the variable is a bit of a misnomer
                     SphAngles(n_pol_out, sphangles);
                     b1.Fresnel_Pol_D[0] = n_pol_out[0];
                     b1.Fresnel_Pol_D[1] = n_pol_out[1];
                     b1.Fresnel_Pol_D[2] = n_pol_out[2];
                     b1.Fresnel_Pol_D_phi = sphangles[0] * RAD2DEG;
                     b1.Fresnel_Pol_D_theta = sphangles[1] * RAD2DEG;
//cout<<"phi theta pol_out deg "<<sphangles[0] * RAD2DEG<<" "<<sphangles[1] * RAD2DEG<<endl;

//n_pol is used further down for the trees, so I'll just equate the correct Fresnel'd n_pol_out to n_pol for now:
                     n_pol[0] = n_pol_out[0];
                     n_pol[1] = n_pol_out[1];
                     n_pol[2] = n_pol_out[2];
                  }//end of if(<400m)

                  nsignal_atAT[0] = sin(theta2) * cos(phi_nposnu2AT);
                  nsignal_atAT[1] = sin(theta2) * sin(phi_nposnu2AT);
                  nsignal_atAT[2] = cos(theta2);

//I think it's ok to use same 'nnu' here, because imagine neutrino interacting in firn, that should give same n_pol
// NO, ABOVE STATEMENT IS WRONG. 08/30/11: KD & SB FOUND OUT THAT IT DOESN'T HOLD EXCEPT FOR COUPLE OF SCENARIOS!!!!!!!
                  //GetPolarization(nnu,nsignal_atAT,n_pol, n_Bfield);


               }//END of if(FIRN)//NOTE: here we have overwritten the Polarization if interaction was in ice, if it was in firn, it stayed the same


               if (PLANEWAVE) {
//calculating distance and time from AT to STcenter, using distance of pt from a plane, where plane is defined by its normal
// =|ax2 + by2 + cy2 +d| / sqrt(a^2+b^2+c^2)
                  double d_AT2ST_denominator = (nsignal_atST[0] * ATCoordinates8[WhichAntenna][0] + nsignal_atST[1] * ATCoordinates8[WhichAntenna][1] +
                                                nsignal_atST[2] * ATCoordinates8[WhichAntenna][2] + (-nsignal_atST[0] * ATCoordinate[0] -
                                                      -nsignal_atST[1] * ATCoordinate[1] - nsignal_atST[2] * ATCoordinate[2]));
                  double d_AT2ST_numerator = sqrt(Square(nsignal_atST[0]) + Square(nsignal_atST[1]) + Square(nsignal_atST[2]));
                  double d_AT2ST = d_AT2ST_denominator / d_AT2ST_numerator ; //why do I have it inverted??
//double time_AT2ST=0; //time offset compared to center of ST

                  double time_AT2ST = d_AT2ST * NICE / C;
                  if (FIRN)
                     time_AT2ST = d_AT2ST * NFIRN / C; //station is always in firn, right?


                  b1.my_time_AT2ST[WhichAntenna] = time_AT2ST;
                  b1.nsignal_i_atAT[WhichAntenna] =  nsignal_atAT[0];
                  b1.nsignal_j_atAT[WhichAntenna] =  nsignal_atAT[1];
                  b1.nsignal_k_atAT[WhichAntenna] =  nsignal_atAT[2];
                  b1.nposnu2AT_i[WhichAntenna] =  nposnu2AT[0]; // CAREFUL ABOUT USAGE IN CASE OF FIRN
                  b1.nposnu2AT_j[WhichAntenna] =  nposnu2AT[1]; // CAREFUL ABOUT USAGE IN CASE OF FIRN
                  b1.nposnu2AT_k[WhichAntenna] =  nposnu2AT[2]; // CAREFUL ABOUT USAGE IN CASE OF FIRN
               }//END of if(PLANEWAVE)


               b1.Bfield_D_firn[0] = n_Bfield[0];
               b1.Bfield_D_firn[1] = n_Bfield[1];
               b1.Bfield_D_firn[2] = n_Bfield[2];
               b1.Polarization_D[0] = n_pol[0];
               b1.Polarization_D[1] = n_pol[1];
               b1.Polarization_D[2] = n_pol[2];


               SphAngles(n_pol, sphangles);
               b1.polariz_phi_D = sphangles[0] * RAD2DEG;
               b1.polariz_theta_D = sphangles[1] * RAD2DEG;


               /*
               if (WhichAntenna==7){
               cout<<
               "inu: "<<inu<<" nnu:"<< nnu[0] <<","<<nnu[1]<<","<<nnu[2]<<
               "n_pol[]=["<<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<"]"<<
               " pol_phi "<<sphangles[0]*RAD2DEG<<"deg, pol_theta"<<sphangles[1]*RAD2DEG<<"deg"<<endl;
               cout<<"--------------"<<endl;
               }
               */

               double sum_vmmhz = 0.;
               double sum_vmmhz1m_unattened = 0; //but including the deviation from the Cherenkov axis
               double sum_vmmhz1m_unattened_max = 0.; //on Cherenkov axis
               double volt_max = 0.;

               Zero(vmmhz, NFREQ);
               Zero(vmmhz1m, NFREQ);
               Zero(vmmhz1m_unattened, NFREQ);

               double   Esumvmmhz1m_unattened = 0.;
               for (int i = 0; i < NFREQ; i++) {
                  vmmhz1m_unattened[i] = GetVmMHz1m(pnu, freq[i], X0ICE, ECICE, n1, AEX_ICE);
//-------------------------------
                  b7.Evmmhz1m_unattened[i] = vmmhz1m_unattened[i];
//-------------------------------
                  vmmhz1m_unattened_max[i] = vmmhz1m_unattened[i]; //KD: are we just taking the one corresponding to max freq? No, we are creating an array that gets updated with the Cerenkov on-cone max further down

                  Esumvmmhz1m_unattened += FREQ_BIN * vmmhz1m_unattened[i];
               }
//-------------------------------
               b7.Esumvmmhz1m_unattened = Esumvmmhz1m_unattened;
//-------------------------------

               // clock_t final_time=clock();
               // cout<<final_time<<endl;


               if (SECKEL == 0) {
                  vmmhz_max = VmMHz(vmmhz1m_max, d_posnu2AT);

                  if (ATTEN_FREQ) { //if attenuation length is a function of frequency
                     for (int i = 0; i < NFREQ; i++) {
                        if (freq[i] <= 122.) {
                           attenlength_up_freq = attenlength_up * (700. - 250.*(freq[i] - 50.) / 72.) / 450.; //freq in MHz, attenlength in m
                           vmmhz_max_freq[i] = VmMHz_attenuated(d_posnu2AT, vmmhz_max, attenlength_up_freq);
                        } else
                           vmmhz_max_freq[i] = VmMHz_attenuated(d_posnu2AT, vmmhz_max, attenlength_up);
                     }
                  }

                  vmmhz_max = VmMHz_attenuated(d_posnu2AT, vmmhz_max, attenlength_up);

                  if (ATTEN_FREQ) {
//                     if (ST_TYPE == 4) {
                        if (vmmhz_max_freq[0]*heff_max_LPA * BW < NV)
                           continue;
//                     }
                     GetVmMHz_freq(vmmhz_max_freq, vmmhz1m_max, pnu, vmmhz);
                  }

                  else {
//                     if (ST_TYPE == 4) {
                        if (vmmhz_max * heff_max_LPA * BW < NV)
                           continue;
//                     }
                     //    cout<<"I'm here 2"<<endl;
                     GetVmMHz(vmmhz_max, vmmhz1m_max, pnu, vmmhz); //KD there's a loop over freq here, creates vmmhz[i]

                     double Esum_vmmhz = 0.;
                     for (int i = 0; i < NFREQ; i++) {
//-------------------------------
                        b7.Evmmhz[WhichAntenna][i] = vmmhz[i];
//-------------------------------

                        Esum_vmmhz += FREQ_BIN * vmmhz[i];
                     }
//-------------------------------
                     b7.Esum_vmmhz[WhichAntenna] = Esum_vmmhz; // =(attenfactor*Esumvmmhz1m_unattened/dis)
//-------------------------------
                  }
                  // Make a vector of V/m/MHz scaled by 1/r and attenuated.
                  // Calculates Jaime's V/m/MHz at 1 m for each frequency
                  // then multiplies by scale factor vmmhz_max/vmmhz1m_max
               }


               for (int i = 0; i < NFREQ; i++) {

                  if (SECKEL == 0) {
                     TaperVmMHz(viewangle, deltheta_em[i], deltheta_had[i], emfrac, hadfrac, vmmhz[i], vmmhz_em[i], vmmhz_had[i]);
//-------------------------------
                     b7.Evmmhz_taper[WhichAntenna][i] = vmmhz[i];
                     b7.Evmmhz_em[WhichAntenna][i] = vmmhz_em[i];
                     b7.Evmmhz_had[WhichAntenna][i] = vmmhz_had[i];
//-------------------------------
                     TaperVmMHz(viewangle, deltheta_em[i], deltheta_had[i], emfrac, hadfrac, vmmhz1m_unattened[i], vmmhz_em[i], vmmhz_had[i]);
//b1.Efield_atposnu_my[WhichAntenna][i] = vmmhz1m_unattened[i];
//-------------------------------
                     b7.Efield_atposnu_my[WhichAntenna][i] = vmmhz1m_unattened[i];
//-------------------------------
                     TaperVmMHz(changle, deltheta_em[i], deltheta_had[i], emfrac, hadfrac, vmmhz1m_unattened_max[i], vmmhz_em[i], vmmhz_had[i]);
                  }

                  if (SECKEL == 1) {
                     vmmhz1m[i] = VmMHz1m(viewangle, freq[i], emfrac, hadfrac, em_shower_length, had_shower_length); //voltage,V/MHz
                     vmmhz[i] = VmMHz(vmmhz1m[i], d_posnu2AT);
                     vmmhz[i] = VmMHz_attenuated(d_posnu2AT, vmmhz[i], attenlength_up); //get vmmhz attenuated
                  }

                  //KD for each of three cases above, sum up. used for speeding up computation
                  sum_vmmhz += FREQ_BIN * vmmhz[i]; //V/m
                  sum_vmmhz1m_unattened += FREQ_BIN * vmmhz1m_unattened[i];
                  sum_vmmhz1m_unattened_max += FREQ_BIN * vmmhz1m_unattened_max[i];

               }//end of freq for loop
//-------------------------------
               b7.Esum_vmmhz_taper[WhichAntenna] = sum_vmmhz;
//-------------------------------

//cout<<"|"<<inu<<"(DIR):"<<iRow_oncone.at(WhichStation)<<","<<iCol_oncone.at(WhichStation)<<","<<NV;

//               if (ST_TYPE == 4) {
                  volt_max = sum_vmmhz * heff_max_LPA;
//cout<<"KD8: "<< "volt_max:" << volt_max << " sum_vmmhz:" << sum_vmmhz << " heff_max_LPA:" << heff_max_LPA << endl;
                  if (volt_max < NV) {
                     b1.iAT[WhichAntenna]    = 0;
                     //b1.kviewangle[WhichAntenna] = 0;
                     //b1.volt_max[WhichAntenna]   = 0;
                     b1.volt_LPA[WhichAntenna]  = 0;
                     //b1.sum_vmmhz[WhichAntenna]  = 0;
                     //b1.sum_vmmhz1m_unattened[WhichAntenna]  = 0;
                     b1.hitangle_e_LPA[WhichAntenna]  = 0;
                     b1.hitangle_h_LPA[WhichAntenna]  = 0;
                     b1.totaltime[WhichAntenna] = 0;
                     continue;
                  }
//               }

               abs_time = (hy1 * NICE + hy2 * NFIRN) / C;

//Replace with new antenna model               if (ST_TYPE == 4) {
/* //OldShelfMC                 
		  if (FIRN)
                     GetHitAngle_LPA(WhichAntenna, N_Ant_perST, nsignal_atAT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA); //KD:added N_Ant_perST
                  else
                     GetHitAngle_LPA(WhichAntenna, N_Ant_perST, nposnu2AT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA);
*/
	       //antenna orientation vectors to pass into Heff function
	          double n_boresight[3];
		  double n_eplane[3];
		  for (int i =0; i<3; i++){
		    n_boresight[i]=Ant_n_boresight[WhichAntenna][i];
		    n_eplane[i]=Ant_n_eplane[WhichAntenna][i];
		  }

		  //Give an e and h plane component for output tree
		  if (FIRN)
		    GetHitAngle(n_boresight, n_eplane, nsignal_atAT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA);
                  else
		    GetHitAngle(n_boresight, n_eplane, nposnu2AT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA);
		  
                  double volt_LPA = 0; //volts of log periodic antenna
                  double volt_LPA_preNoise = 0;
                  term_LPA = 0; //zero term_LPA
                  term_LPA_e = 0;
                  term_LPA_h = 0;

                  b1.e_component_LPA[WhichAntenna] = e_component_LPA;
                  b1.h_component_LPA[WhichAntenna] = h_component_LPA;

                  for (int i = 0; i < NFREQ; i++) { //here needs to be modified

//term_LPA=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(gainv,freq[i]*1.E6)*
//sqrt(pow(e_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[0][i])*(hitangle_e_LPA/flare[0][i]))*exp(-2*ALOG2*(hitangle_h_LPA/flare[1][i])*(hitangle_h_LPA/flare[1][i])),2));

//06/20/11: SB: an average as to get G=G'G(0) where G'=0.5(G_e + G_h)

		    /* //2-13-2017 Moved to GetHeff()
                     term_LPA = vmmhz[i] * FREQ_BIN * 0.5 * GaintoHeight(gainv, freq[i] * 1.E6) *
                                sqrt((pow(e_component_LPA * exp(-2 * ALOG2 * (hitangle_e_LPA / flare[0][i]) * (hitangle_e_LPA / flare[0][i])), 2)  +   pow(e_component_LPA * exp(-2 * ALOG2 * (hitangle_h_LPA / flare[1][i]) * (hitangle_h_LPA / flare[1][i])), 2)) / 2);
		    */
		    
		    if (FIRN)
		      term_LPA = vmmhz[i] * FREQ_BIN * 0.5 * GetHeff(AntType[WhichAntenna], freq[i],n_boresight,n_eplane, nsignal_atAT, n_pol);
		    else
		      term_LPA = vmmhz[i] * FREQ_BIN * 0.5 * GetHeff(AntType[WhichAntenna], freq[i],n_boresight,n_eplane, nposnu2AT, n_pol);


                     /*//works for gain tests
                     term_LPA=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(gainv,freq[i]*1.E6)*sqrt(pow(e_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[0][i])*(hitangle_e_LPA/flare[0][i])),2)     +
                     0.01*pow(h_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[1][i])*(hitangle_e_LPA/flare[1][i])),2));
                     */


                     /*
                     //KD: added 6/9/11 for H-plane averaging
                     //is gainv the same?
                     //term_LPA_h=term_LPA_e; 06/13/11 I had equation below but we want to average the exponetional term
                     term_LPA_h=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(gainv,freq[i]*1.E6)*sqrt(pow(h_component_LPA*exp(-2*ALOG2*(hitangle_h_LPA/flare[1][i])*(hitangle_h_LPA/flare[1][i])),2) + 0.01*pow(e_component_LPA*exp(-2*ALOG2*(hitangle_h_LPA/flare[0][i])*(hitangle_h_LPA/flare[0][i])),2));
                     //DO THE AVERAGE:
                     term_LPA = 0.5*(term_LPA_e + term_LPA_h);

                     */

                     /*
                     //old gain EQUATION
                     term_LPA=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(GetGainV(freq[i]*1.E6),freq[i]*1.E6)*sqrt(pow(e_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[0][i])*(hitangle_e_LPA/flare[0][i])),2)
                                                 +
                     0.01*pow(h_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[1][i])*(hitangle_e_LPA/flare[1][i])),2));
                     */

                     //-----------------------
                     b1.heff_my[i]  = GaintoHeight(gainv, freq[i] * 1.E6); //KD 02/03/2011 I have taken out antenna dependence, as there's none //b1.heff_my[WhichAntenna][i] = GaintoHeight((GetGainV(freq[i])),freq[i]*1.E6);
                     b1.vmmhz_my[WhichAntenna][i]  = vmmhz[i];
                     b1.term_LPA_my[WhichAntenna][i]  = term_LPA;

//-------------------------------
                     b7.heff_my[i]  = GaintoHeight(gainv, freq[i] * 1.E6);
                     b7.vmmhz_my[WhichAntenna][i]  = vmmhz[i];
                     b7.term_LPA_my[WhichAntenna][i]  = term_LPA;
//-------------------------------


//cout<<"i "<<i<<" freq "<< freq[i]<< " term_LPA[i] "<< term_LPA <<" vmmhz[i] " << vmmhz[i] << endl;
                     volt_LPA += term_LPA;

//-------------
//cout<<" "<<gain<<"  "<<endl;
//cout<<"i="<<i<<" gainv_measured= "<<gainv_measured[i]<<endl;

//cout<<"i="<<i<<" gainv = "<<gainv<<" flares ="<< flare[0][i] <<" , "<< flare[1][i] <<" , "<< flare[2][i]<<" , "<< flare[3][i] <<endl;

//cout<<"i:"<< i <<",frequency:"<< (freq[i]*1.E6) <<",GetGainV:" << GetGainV(freq[i]*1.E6) <<endl;
//cout<<"flares ="<< flare[0][i] <<" , "<< flare[1][i] <<endl;
//cout<<"posnu"<< posnu[1] <<","<< posnu[2] <<","<<posnu[3]<<endl;
//-------------

                  }
//cout<<"--------------------------"<<endl;

//if (eshower_em > 1.E10){
//cout<<"inu= "<<inu<<" elpm= "<< elpm <<" vmmhz1m_max= "<<vmmhz1m_max<<endl;
//cout<<"inu= "<<inu<<" pnu, em-, had-shower= ("<<pnu<<","<<eshower_em<<","<<eshower_had<<") volt_max, -LPA= ("<<volt_max<<","<<volt_LPA<<") changle, viewangle= ("<<changle*RAD2DEG<<","<<viewangle*RAD2DEG<<") sum_vmmhz ,-unattened, -unattened_max= ("<<sum_vmmhz<<","<<sum_vmmhz1m_unattened<<","<<sum_vmmhz1m_unattened_max<<")"<<endl;
//}

//outantposnu<< posnu[0] <<" "<< posnu[1] <<" "<<posnu[2]<<" "<<nnu[0] <<" "<< nnu[1] <<" "<<nnu[2]<< " " << weight <<" "<< ATCoordinate[0]<<" "<<ATCoordinate[1]<<" "<<ATCoordinate[2]<< endl;

//introduced 08/01/11


                  if (SIGNAL_FLUCT) {  // KD: 10/11/10 adding time fluctation here
                     volt_LPA_preNoise = volt_LPA; //gets populated if Noise present
                     volt_LPA += Rand3.Gaus(0., VNOISE);
//cout<<inu<<" " <<Rand3.Gaus(0.,VNOISE)<<endl;
                     abs_time += Rand3.Gaus(0., TIMEPRECISION);
                  }

                  if (SHADOWING) {
                     if (FIRN) {
                        //KD adding shadowing cut only for direct
                        if (posnu[2] > (ICETHICK - FIRNDEPTH)) { //here means it's in firn

                           if (sqrt((posnu[0] - ATCoordinate[0]) * (posnu[0] - ATCoordinate[0]) + (posnu[1] - ATCoordinate[1]) * (posnu[1] - ATCoordinate[1])) > (GetRange(posnu[2]) + 20.4)) //adding 25.7 or 20.4 to allow for further bending up to -2m. Also making it with ref to station center coord, for multi station simulation.
                              continue;
                        }

                        if (posnu[2] < (ICETHICK - FIRNDEPTH)) { //here in ice layer
//making ice-firn boundary ray more dynamic 11/27/11
                           if (sqrt((posnu[0] - ATCoordinate[0]) * (posnu[0] - ATCoordinate[0]) + (posnu[1] - ATCoordinate[1]) * (posnu[1] - ATCoordinate[1])) > (GetRange(ICETHICK - FIRNDEPTH) + 20.4 + (ICETHICK - FIRNDEPTH - posnu[2])))

                              //if (sqrt((posnu[0]- ATCoordinate[0])*(posnu[0]- ATCoordinate[0]) + (posnu[1]- ATCoordinate[1])*(posnu[1]- ATCoordinate[1])) > (175.9 + 20.4 + (ICETHICK-FIRNDEPTH-posnu[2])))
                              // we took radial distance as [firn distance + the distance below firn] as we approximated 45degrees for shadowing, as this is nearly the TIR angle from 1.8 to 1.3
                              continue;
                        }
                     }
                  }

//cout<<"||"<<inu<<"(DIR):"<<iRow_oncone.at(WhichStation)<<","<<iCol_oncone.at(WhichStation)<<","<<NV<<endl;

                  if (volt_LPA >= NV) {
                     //cout<<"I'm here 6.0"<<endl;
                     sum_triggeredAT++;
                     iAT_L1.push_back(WhichAntenna);
                     path_inice_L1.push_back(hy1);
                     path_infirn_L1.push_back(hy2);
                     totaltime_L1.push_back(abs_time);
                     Efield_atposnu_L1.push_back(sum_vmmhz1m_unattened);
                     volte_allAT_L1.push_back(volt_LPA);
                     volth_allAT_L1.push_back(0);
                     AT_posi_x_L1.push_back(ATCoordinates8[WhichAntenna][0]);
                     AT_posi_y_L1.push_back(ATCoordinates8[WhichAntenna][1]);
                     AT_posi_z_L1.push_back(ATCoordinates8[WhichAntenna][2]);
                     Dir_shower_x_L1.push_back(nsignal_atAT[0]);
                     Dir_shower_y_L1.push_back(nsignal_atAT[1]);
                     Dir_shower_z_L1.push_back(nsignal_atAT[2]);
                     theta_shower_atposnu_L1.push_back(theta1);
                     //    cout<<"I'm here 4.6"<<endl;
                     hitangle_e_LPA_L1.push_back(hitangle_e_LPA * RAD2DEG);
                     hitangle_h_LPA_L1.push_back(hitangle_h_LPA * RAD2DEG);

//b1.kviewangle[WhichAntenna]    = viewangle*RAD2DEG;
//b1.volt_max[WhichAntenna]      = volt_max;
//b1.volt_LPA[WhichAntenna]      = volt_LPA;
//b1.sum_vmmhz[WhichAntenna]     = sum_vmmhz;
//b1.sum_vmmhz1m_unattened[WhichAntenna]  = sum_vmmhz1m_unattened;
//b1.hitangle_e_LPA[WhichAntenna]      = hitangle_e_LPA*RAD2DEG;

                  }

//cout<<"KD10: "<<"WhichAntenna:"<<WhichAntenna<<" volt_LPA=" << volt_LPA << "  NV=" << NV << endl;
//cout<<endl;

//b1.kviewangle[WhichAntenna]    = viewangle*RAD2DEG;
//b1.volt_max[WhichAntenna]      = volt_max;
                  b1.volt_LPA[WhichAntenna]     = volt_LPA;
                  b1.volt_LPA_preNoise[WhichAntenna]  = volt_LPA_preNoise;
//b1.sum_vmmhz[WhichAntenna]     = sum_vmmhz;
//b1.sum_vmmhz1m_unattened[WhichAntenna]  = sum_vmmhz1m_unattened;
//b1.hitangle_e_LPA[WhichAntenna]      = hitangle_e_LPA*RAD2DEG;

                  //  cout<<"I'm here 4.7"<<endl;
//replace with new antenna model               }
               b1.sum_vmmhz1m_unattened_max  = sum_vmmhz1m_unattened_max;
//b1.Emax_atposnu       = sum_vmmhz1m_unattened_max;
               b1.hitangle_e_LPA[WhichAntenna]     = hitangle_e_LPA * RAD2DEG;
               b1.hitangle_h_LPA[WhichAntenna]     = hitangle_h_LPA * RAD2DEG;
               b1.theta_my_signalAT[WhichAntenna]  = theta2 * RAD2DEG;
               b1.phi_my_signalAT[WhichAntenna] = phi_nposnu2AT * RAD2DEG;
               b1.totaltime[WhichAntenna]    = abs_time;
//b1.volt_LPA[WhichAntenna]      = volte_allAT_L1.at(volte_allAT_L1.size()-1);
//cout<<" KD2: "<< inu<<" , " <<abs_time <<endl ;       //  cout<<"I'm here 4.8"<<endl;
            }//end of for antennas per station loop

//cout<<"=======================" << endl;

//++cout <<"KD  : Having populated L1 arrays, I can now output L1 details:"<<endl;

//++for (int i=0; i<iAT_L1.size(); i++) //KD: looping over triggered antenna
//++{
//++cout<<"KD11a: "<<"iAT_L1[i]:"<< iAT_L1[i] << "  AT position L1[i]=("<< AT_posi_x_L1[i] << ","<< AT_posi_y_L1[i] << "," << AT_posi_z_L1[i] <<")  volt_LPA="<< volte_allAT_L1[i]<<" totaltime_L1="<<totaltime_L1[i]<<endl;
//"  Efield_atposnu_L1[i]="<< Efield_atposnu_L1[i] << endl;
//cout<<"KD9b: "<< "AT_posi_x_L1[i]="<< AT_posi_x_L1[i] << "  AT_posi_y_L1[i]="<< AT_posi_y_L1[i] << "  AT_posi_z_L1[i]=" << AT_posi_z_L1[i] << endl;
//++}
//cout<<"min time L1"<<Min(totaltime_L1,iAT_L1.size())<<endl; //KD Trying a smart to get Min and Max for time differences

//++cout<<"KD11b: sum_triggeredAT: " << sum_triggeredAT << endl;
//cout<<"------------------------------------------------------------------------------"<<endl;




//02/2011: above we have evaluated all antennas and will now consider trigger scenarios

            if (sum_triggeredAT >= N_Ant_Trigger - 1) //KD : new variable from input file
//       if(sum_triggeredAT>=2)    //if the station is triggered(L2 satisfied)
            {
               // START OF TRIGGER SATISFIED ANTENNAS CONDITION

               // cout<<sum_triggeredAT<<endl;
//++cout<<"KD  : I'm here in a loop since sum_triggeredAT>=2"<<endl;

//KD 02/17/11: THE ABOVE IS AN EXTRA STEP AS I EVENTUALLY DISCARD THE '==n-1' CONDITION IMMEDIATELY BELOW
//               if (ST_TYPE == 4) {
                  if (sum_triggeredAT == N_Ant_Trigger - 1) //KD: new variable from input file
//          if(sum_triggeredAT==2)//at least 3 out of 8 LPAs
                     continue;   //KD: it looks like here, I just terminate surface stations calculations and move on to mirror calculation.
//               }
//++cout<<"KD  : I'm continuing presumably because sum_t'AT isn't 2 but is: "<<sum_triggeredAT <<endl;


               if (0) { //activate for studying different scenario studies
//-------------------------------------------------------------
//-----STUART's question---------------------------------------
                  int even_ant = 0; //resetting variable
                  int odd_ant = 0; //resetting variable

                  for (uint WhichAnt = 0; WhichAnt < iAT_L1.size(); WhichAnt++) { //go through all triggered AT sequentially
                     //if(iAT_L1[WhichAnt]==0 || iAT_L1[WhichAnt]==2 || iAT_L1[WhichAnt]==4 ||iAT_L1[WhichAnt]==6)
                     if (iAT_L1[WhichAnt] == 0) { //ignoring AT=2   || iAT_L1[WhichAnt]==4
                        even_ant++;//incrementing count for even antennas
                     }

                     if (iAT_L1[WhichAnt] == 1 || iAT_L1[WhichAnt] == 3 || iAT_L1[WhichAnt] == 5) { //ignoring
                        odd_ant++;//incrementing count for even antennas
                     }
                  }
//after going through all loops and getting even and odd triggered number of ant's
                  if (even_ant < 1 && odd_ant < 2) {
                     continue; //discard if not 1 AND 3, go out of loop and do not trigger station
                  }
//------------------------------------------------------------
               } //end of different scenario


               //  cout<<"I'm here 5"<<endl;
               sum_triggeredST++;
               N_TriggeredAT_perST_L2.push_back(sum_triggeredAT);
               Trig_Type_L2.push_back(1);
               iRow_L2.push_back(iRow_oncone.at(WhichStation));
               iCol_L2.push_back(iCol_oncone.at(WhichStation));
//++cout<<"KD  : The code has just updated 'sum_triggeredST' and filled first four L2 vectors."<<endl;


               for (int i = 0; i < sum_triggeredAT; i++) {
                  iAT_L2.push_back(iAT_L1.at(i));
                  path_inice_L2.push_back(path_inice_L1.at(i));
                  path_infirn_L2.push_back(path_infirn_L1.at(i));
                  totaltime_L2.push_back(totaltime_L1.at(i));
                  Efield_atposnu_L2.push_back(Efield_atposnu_L1.at(i));
                  volte_allAT_L2.push_back(volte_allAT_L1.at(i));
                  volth_allAT_L2.push_back(volth_allAT_L1.at(i));
                  AT_posi_x_L2.push_back(AT_posi_x_L1.at(i));
                  AT_posi_y_L2.push_back(AT_posi_y_L1.at(i));
                  AT_posi_z_L2.push_back(AT_posi_z_L1.at(i));
                  Dir_shower_x_L2.push_back(Dir_shower_x_L1.at(i));
                  Dir_shower_y_L2.push_back(Dir_shower_y_L1.at(i));
                  Dir_shower_z_L2.push_back(Dir_shower_z_L1.at(i));
                  theta_shower_atposnu_L2.push_back(theta_shower_atposnu_L1.at(i));
//++cout<<"KD  : This repeats sum_t'AT times as I fill up remaining fourteen L2 vectors.+"<<endl;
                  hitangle_e_LPA_L2.push_back(hitangle_e_LPA_L1.at(i));
                  hitangle_h_LPA_L2.push_back(hitangle_h_LPA_L1.at(i));
               }

//KD: Should be filled here within a loop as it passes AT/ST test
               double Efield_atposnu_L2_mymax   = 0;
               double volte_allAT_L2_mymax   = 0;
               double volte_allAT_L2_mymin = 31.55; //arbitrary minimum
               double hitangle_e_LPA_mymax   = 0;
               double hitangle_h_LPA_mymax   = 0;
               double totaltime_L2_mymax  = 0;
               double totaltime_L2_mymin  = 31.55; //arbitrary max

               for (uint j = 0; j < iAT_L2.size(); j++) { //KD: 10/9/10 only within triggered data
//cout<<"iAT = "<< iAT_L2.at(j)<< " Efield_atposnu= " << Efield_atposnu_L2.at(j) << " volt_e= "<< volte_allAT_L1.at(j)<<endl;
                  b1.iAT[j]   = iAT_L2[j] + 1 ;

                  if (Efield_atposnu_L2[j] > Efield_atposnu_L2_mymax)
                     Efield_atposnu_L2_mymax = Efield_atposnu_L2[j];

                  if (volte_allAT_L2[j] > volte_allAT_L2_mymax) {
                     volte_allAT_L2_mymax = volte_allAT_L2[j];
                     hitangle_e_LPA_mymax = hitangle_e_LPA_L2[j];
                     hitangle_h_LPA_mymax = hitangle_h_LPA_L2[j];
                     iAT_index_mymax = iAT_L2[j];
                  }

                  //selecting minimum voltage and its index
                  if (volte_allAT_L2[j] < volte_allAT_L2_mymin) {
                     volte_allAT_L2_mymin = volte_allAT_L2[j];
                     iAT_index_mymin = iAT_L2[j];
                  }

                  if (totaltime_L2[j] > totaltime_L2_mymax) {
                     totaltime_L2_mymax = totaltime_L2[j];
                     iAT_index_timemax = iAT_L2[j];
                  }

                  if (totaltime_L2[j] < totaltime_L2_mymin) {
                     totaltime_L2_mymin = totaltime_L2[j];
                     iAT_index_timemin = iAT_L2[j];
                  }
               }

//Efield_atposnu_L2_mymax = Max(Efield_atposnu_L2, Efield_atposnu_L2.size());
//cout<<"Efield_atposnu_L2_mymax= "<< Efield_atposnu_L2_mymax <<" volte_allAT_L2_mymax= "<<volte_allAT_L2_mymax<<endl;
               b1.Efield_atposnu_mymax    = Efield_atposnu_L2_mymax;
               b1.volte_allAT_mymax    = volte_allAT_L2_mymax;
               b1.volte_allAT_mymin        = volte_allAT_L2_mymin;
               b1.hitangle_e_LPA_mymax    = hitangle_e_LPA_mymax;
               b1.hitangle_h_LPA_mymax    = hitangle_h_LPA_mymax;
               b1.iAT_index_mymax      = iAT_index_mymax;
               b1.iAT_index_mymin      = iAT_index_mymin;
               b1.totaltime_mymin      = totaltime_L2_mymin;
               b1.totaltime_mymax      = totaltime_L2_mymax;
               b1.iAT_index_timemax    = iAT_index_timemax;
               b1.iAT_index_timemin    = iAT_index_timemin;
//b1.viewangle[WhichAntenna]     = viewangle*RAD2DEG;
//b1.volt_max[WhichAntenna]      = volt_max;
//b1.volt_LPA[WhichAntenna]      = volt_LPA;
//b1.sum_vmmhz[WhichAntenna]     = sum_vmmhz;
//b1.sum_vmmhz1m_unattened[WhichAntenna]  = sum_vmmhz1m_unattened;
            }  //END OF TRIGGER SATISFIED ANTENNAS CONDITION




//b1.sum_vmmhz1m_unattened_max   = sum_vmmhz1m_unattened_max;

//b1.Emax_atposnu       = sum_vmmhz1m_unattened_max;
//b1.E_atposnu       = Efield_atposnu_L1;    //KD: should be each AT
//b1.Abs_time        = totaltime_L1;   //KD: should be each AT
//b1.Dir_Shower         = theta_shower_atposnu_L1; //KD: should be each AT
//b1.Path_inice         = path_inice_L1;  //KD: should be each AT
//b1.Path_infirn        = path_infirn_L1; //KD: should be each AT



//++cout <<"KD  : I'm now outside the 'sum_AT>=2' loop with sum_triggeredAT still: "<<sum_triggeredAT<<endl;
//++cout <<"KD  : I can now output L2 details, that should be exact as L1 above:"<<endl;
//++for (int i=0; i<sum_triggeredAT; i++) //KD: looping over triggered antenna, was wrongly using iAT_L2.size previously??
//++{
//++cout<<"KD12a: "<<"iAT_L2[i]:"<< iAT_L2[i] << "  AT position L2[i]=("<< AT_posi_x_L2[i] << ","<< AT_posi_y_L2[i] << "," << AT_posi_z_L2[i] <<")  volt_LPA="<< volte_allAT_L2[i]<< " totaltime_L2="<<totaltime_L1[i]<<endl;
//"  Efield_atposnu_L2[i]="<< Efield_atposnu_L2[i] << endl;
//++}
//++cout<<"KD12b: sum_triggeredAT:"<< sum_triggeredAT <<endl;
//++cout<<"KD12c: sum_triggeredST:"<< sum_triggeredST <<endl;

         }//end of surface stations loop

//++cout <<"KD  : HERE I'M NOW MOVING ON TO MIRROR CALCULATIONS"<<endl;
//==============================================================================
//--------------------MIRROR----------------------------------------------------
//==============================================================================

         // cout<<"I'm here 6"<<endl;
         //start the mirror stations
         //variables after station trigger(L2)
         vector<int> N_TriggeredAT_perST_mirror_L2;//the number of triggered antennas in each triggered station
         vector<int> Trig_Type_mirror_L2;//triggered by direct or reflected signal;
         vector<int> iRow_mirror_L2;
         vector<int> iCol_mirror_L2;

         //into each triggered station, analyze the antennas in the station
         vector<int> iAT_mirror_L2;//in the triggered station, which antennas are triggered, this is the label of the triggered antenna
         vector<double> path_inice_mirror_L2;
         vector<double> path_infirn_mirror_L2;
         vector<double> totaltime_mirror_L2;//the time for a signal getting to the antenna
         vector<double> Efield_atposnu_mirror_L2;
         vector<double> volte_allAT_mirror_L2;
         vector<double> volth_allAT_mirror_L2;
         vector<double> AT_posi_x_mirror_L2;
         vector<double> AT_posi_y_mirror_L2;
         vector<double> AT_posi_z_mirror_L2;
         vector<double> Dir_shower_x_mirror_L2;
         vector<double> Dir_shower_y_mirror_L2;
         vector<double> Dir_shower_z_mirror_L2;
         vector<double> theta_shower_atposnu_mirror_L2;
         vector<double> hitangle_e_LPA_mirror_L2;
         vector<double> hitangle_h_LPA_mirror_L2;

         N_TriggeredAT_perST_mirror_L2.clear();
         Trig_Type_mirror_L2.clear();
         iRow_mirror_L2.clear();
         iCol_mirror_L2.clear();


         iAT_mirror_L2.clear();
         path_inice_mirror_L2.clear();
         path_infirn_mirror_L2.clear();
         totaltime_mirror_L2.clear();
         Efield_atposnu_mirror_L2.clear();
         volte_allAT_mirror_L2.clear();
         volth_allAT_mirror_L2.clear();
         AT_posi_x_mirror_L2.clear();
         AT_posi_y_mirror_L2.clear();
         AT_posi_z_mirror_L2.clear();
         Dir_shower_x_mirror_L2.clear();
         Dir_shower_y_mirror_L2.clear();
         Dir_shower_z_mirror_L2.clear();
         theta_shower_atposnu_mirror_L2.clear();
         hitangle_e_LPA_mirror_L2.clear();
         hitangle_h_LPA_mirror_L2.clear();

//----KD added 7/8
         //for plotting
         viewangle_triggered.clear();
         viewangle_triggered_mirror.clear();
         attenfactor.clear();
         attenfactor_mirror.clear();
         dis.clear();
         dis_mirror.clear();
         vdis.clear();
         vdis_mirror.clear();
         cosz.clear();
         cosz_mirror.clear();
         coszp.clear();
         coszp_mirror.clear();

         theta_signal.clear();
         theta_signal_atAT.clear();
         theta_signal_mirror.clear();
         theta_signal_atAT_mirror.clear();
//---from ST<3--------



         sum_triggeredST_mirror = 0;


//KD start of mirror stations loop for st_type=3 or 4
         for (uint WhichMirrorStation = 0; WhichMirrorStation < iRow_oncone_mirror.size(); WhichMirrorStation++) {

            //cout<<"4th place"<<double(clock()-init_time)/(double)CLOCKS_PER_SEC<<endl;
            //here the GetMIrrorATLocation actually is to get the station location
            GetMirrorATLocation(iRow_oncone_mirror.at(WhichMirrorStation), iCol_oncone_mirror.at(WhichMirrorStation), MirrorATCoordinate);
            // cout<<"line 2779"<<endl;

//-----------------------------------------------------------------------------------------------------------
// ESTABLISHING MAIN TRIGGER + LOOSE TRIGGER CONDITIONS FOR MULTI-STATION
//-----------------------------------------------------------------------------------------------------------
//cout<<inu<<"(M):"<<iRow_oncone_mirror.at(WhichMirrorStation)<<","<<iCol_oncone_mirror.at(WhichMirrorStation)<<endl; // print out to establish varying station row/col
            if (HEXAGONAL) {
               if (iRow_oncone_mirror.at(WhichMirrorStation) == 1 && iCol_oncone_mirror.at(WhichMirrorStation) == 2) {
                  NV = VNOISE * NSIGMA;
               } else {
                  NV = VNOISE * NSIGMA; //Note: JCH changes this from 5.0 to NSIGMA until I can figure out why the threshold is sometimes 5 sigma and sometimes the desired value.
               }
            }
//cout<<inu<<"(REF):"<<iRow_oncone_mirror.at(WhichMirrorStation)<<","<<iCol_oncone_mirror.at(WhichMirrorStation)<<","<<NV;
//-----------------------------------------------------------------------------------------------------------

            //Set Antenna Positions
            double MirrorATCoordinates8[N_Ant_perST][3];//the detailed position of the center of each LPA in a station                        
            if (StationType == 0){ //All antennas pointing down (up), equally spaced around station center
                for (int i = 0; i < N_Ant_perST; i++) {
                   double phi = (2. / N_Ant_perST) * PI * i; //the phi angle of each LPA's center
                   MirrorATCoordinates8[i][0] = MirrorATCoordinate[0] + ST4_R * cos(phi);
                   MirrorATCoordinates8[i][1] = MirrorATCoordinate[1] + ST4_R * sin(phi);
                   MirrorATCoordinates8[i][2] = MirrorATCoordinate[2];
                }
            }
            else if (StationType == 1){ //Some custom antenna config
                cout<<"Haven't defined this yet"<<endl;
            }
            else {
                cout<<"Invalid Station type "<<StationType<<endl;
            }
            
            //Set Antenna Types (This better be the same as for the non-mirrored case)
            int AntType[N_Ant_perST];
            //type 0 = 100MHz theoretical LPDA (original ShelfMC model)
            //type 1 = 100MHz Create LPDA, Anna's WhippleD model
            if (StationType == 0){
                for (int i = 0; i < N_Ant_perST; i++) {
                    AntType[i]=0;
                }
            }
            else if (StationType == 1){
                for (int i = 0; i < N_Ant_perST; i++) {
                    AntType[i]=1;
                }
            }
            else {
                cout<<"Invalid Station type "<<StationType<<endl;
            }

            //Set Antenna Orientation (Flip z component WRT non-mirrored case)
            double MirrorAnt_n_boresight[N_Ant_perST][3];
            double MirrorAnt_n_eplane[N_Ant_perST][3];

            if (StationType == 0) {
                for (int i = 0; i < N_Ant_perST; i++) {
                    MirrorAnt_n_eplane[i][0] = cos((0.5 + i * (2. / N_Ant_perST))*PI);
                    MirrorAnt_n_eplane[i][1] = sin((0.5 + i * (2. / N_Ant_perST))*PI);
                    MirrorAnt_n_eplane[i][2] = 0.;
                    MirrorAnt_n_boresight[i][0] = 0.;
                    MirrorAnt_n_boresight[i][1] = 0.;
                    MirrorAnt_n_boresight[i][2] = 1.; //facing up now, since Stn is mirrored
                }
            }
            else if (StationType == 1) {
                cout<<"haven't defined this yet"<<endl;
            }
            else {
                cout<<"Invalid Station Type"<<endl;
            }

            //define some variables before going into the antenna loop of one station

            //variables after antenna trigger(L1)

            //into each triggered station, analyze the antennas in the station
            vector<int> iAT_mirror_L1;//in the triggered station, which antennas are triggered, this is the label of the triggered antenna
            vector<double> path_inice_mirror_L1;
            vector<double> path_infirn_mirror_L1;
            vector<double> totaltime_mirror_L1;//the time for a signal getting to the antenna
            vector<double> Efield_atposnu_mirror_L1;
            vector<double> volte_allAT_mirror_L1;
            vector<double> volth_allAT_mirror_L1;
            vector<double> AT_posi_x_mirror_L1;
            vector<double> AT_posi_y_mirror_L1;
            vector<double> AT_posi_z_mirror_L1;
            vector<double> Dir_shower_x_mirror_L1;
            vector<double> Dir_shower_y_mirror_L1;
            vector<double> Dir_shower_z_mirror_L1;
            vector<double> theta_shower_atposnu_mirror_L1;
            vector<double> hitangle_e_LPA_mirror_L1;
            vector<double> hitangle_h_LPA_mirror_L1;

            iAT_mirror_L1.clear();
            path_inice_mirror_L1.clear();
            path_infirn_mirror_L1.clear();
            totaltime_mirror_L1.clear();
            Efield_atposnu_mirror_L1.clear();
            volte_allAT_mirror_L1.clear();
            volth_allAT_mirror_L1.clear();
            AT_posi_x_mirror_L1.clear();
            AT_posi_y_mirror_L1.clear();
            AT_posi_z_mirror_L1.clear();
            Dir_shower_x_mirror_L1.clear();
            Dir_shower_y_mirror_L1.clear();
            Dir_shower_z_mirror_L1.clear();
            theta_shower_atposnu_mirror_L1.clear();
            hitangle_e_LPA_mirror_L1.clear();
            hitangle_h_LPA_mirror_L1.clear();



            // cout<<"4.2th place"<<double(clock()-init_time)/CLOCKS_PER_SEC<<endl;
            int sum_triggeredAT_mirror = 0;
            for (int WhichMirrorAntenna = 0; WhichMirrorAntenna < N_Ant_perST; WhichMirrorAntenna++)
//KD  modifying above for 2 out of 3 case
//for(int WhichMirrorAntenna=0; WhichMirrorAntenna<3; WhichMirrorAntenna++)
            {

//               if (ST_TYPE == 4)
               VectorMinus(MirrorATCoordinates8[WhichMirrorAntenna], posnu, new_posnu2MirrorAT);

               nVector(new_posnu2MirrorAT, nposnu2MirrorAT);

               if (nposnu2MirrorAT[1] >= 0) {
                  if (nposnu2MirrorAT[0] > 0)
                     phi_nposnu2MirrorAT = atan(nposnu2MirrorAT[1] / nposnu2MirrorAT[0]);
                  else if (nposnu2MirrorAT[0] < 0)
                     phi_nposnu2MirrorAT = PI + atan(nposnu2MirrorAT[1] / nposnu2MirrorAT[0]);
                  else
                     phi_nposnu2MirrorAT = PI / 2.;
               } else {
                  if (nposnu2MirrorAT[0] > 0)
                     phi_nposnu2MirrorAT = 2 * PI + atan(nposnu2MirrorAT[1] / nposnu2MirrorAT[0]);
                  else if (nposnu2MirrorAT[0] < 0)
                     phi_nposnu2MirrorAT = PI + atan(nposnu2MirrorAT[1] / nposnu2MirrorAT[0]);
                  else
                     phi_nposnu2MirrorAT = 3 * PI / 2.;
               }

               theta_nposnu2MirrorAT = acos(nposnu2MirrorAT[2]); //greater than 90 degrees
               d_posnu2MirrorAT = Mag(new_posnu2MirrorAT);
               viewangle_mirror = Angle(nnu, nposnu2MirrorAT);
               GetPolarization(nnu, nposnu2MirrorAT, n_pol, n_Bfield);
               theta2_mirror = PI - theta_nposnu2MirrorAT; //the angle between nsignal and -z direction, which is less than 90 degrees
               double theta1_mirror = theta2_mirror;
               hy1_mirror = d_posnu2MirrorAT;
               hy2_mirror = 0.; //if no firn
               //  cout<<"here line 2855"<<endl;

               //KD 01/18/11: added these 3 lines for NO FIRN case for PLANEWAVE calcs below
               nsignal_mirror_atAT[0] = sin(theta2_mirror) * cos(phi_nposnu2MirrorAT);
               nsignal_mirror_atAT[1] = sin(theta2_mirror) * sin(phi_nposnu2MirrorAT);
               nsignal_mirror_atAT[2] = -cos(theta2_mirror);

               /*
               cout<<"ievt "<<inu<<" nnu:"<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; LOS P:"
                  <<nsignal_mirror_atAT[0]<<","<<nsignal_mirror_atAT[1]<<","<<nsignal_mirror_atAT[2]<<" angles:"
                  <<theta2_mirror*RAD2DEG<<","<<phi_nposnu2MirrorAT*RAD2DEG<<" >"<<viewangle_mirror*RAD2DEG<<" pol:"
                  <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]
                  <<endl;
               */

               GetBouncePoint(theta_nposnu2MirrorAT, nposnu2MirrorAT, posnu, bounce);
               b1.bouncept[0] = bounce[0];
               b1.bouncept[1] = bounce[1];
               b1.bouncept[2] = bounce[2];

//KD:06/30/11 NEED TO IMPROVE FOR THE SMALL SUBSET OF FIRN EVENTS THAT PROPAGATE THROUGH THE ICE BULK AT SLIGHTLY DIFFERENT ANGLES AND DIRECTIONS


//SHADOWING ONLY OCCURS FOR FIRN PRESENCE
               if (SHADOWING) {
                  if (FIRN) {
                     //KD adding shadowing cut only for reflected

                     /* //HAVE TO THINK HOW TO INSTATE THIS FIRN MODIFIED THETA
                     if(posnu[2]>(ICETHICK-FIRNDEPTH)){ //here means it's in firn
                        if (sqrt((posnu[0] - ATCoordinate[0])*(posnu[0]-ATCoordinate[0]) + (posnu[1]-ATCoordinate[1])*(posnu[1]-ATCoordinate[1])) > (GetRange(posnu[2])+20.4)) //adding 25.7 or 20.4 to allow for further bending up to -2m. Also making it with ref to station center coordinates with ATCoordinate, for multi station simulation.
                        continue;
                              }
                            */


//if(posnu[2]<(ICETHICK-FIRNDEPTH)){ //here in ice layer

                     if (sqrt((bounce[0] - MirrorATCoordinate[0]) * (bounce[0] - MirrorATCoordinate[0]) + (bounce[1] - MirrorATCoordinate[1]) * (bounce[1] - MirrorATCoordinate[1])) > (GetRange(ICETHICK - FIRNDEPTH) + 20.4 + (ICETHICK - FIRNDEPTH - bounce[2])))
//if (sqrt((bounce[0]- MirrorATCoordinate[0])*(bounce[0]- MirrorATCoordinate[0]) + (bounce[1]- MirrorATCoordinate[1])*(bounce[1]- MirrorATCoordinate[1])) > (175.9 + 20.4 + (ICETHICK-FIRNDEPTH-bounce[2])))
// we took radial distance as [firn distance + the distance below firn] as we approximated 45degrees for shadowing, as this is nearly the TIR angle from 1.8 to 1.3. hence horizontal range is same as depth.
                        continue;
                     //      }
                  }//FIRN within SHADOWING
               } // SHADOWING



               if (FIRN) {

                  //NOTE: order of firn calculation is reversed as compared to direct
                  if (posnu[2] < (ICETHICK - FIRNDEPTH)) { //interact in the ice
                     double x1_mirror = 1.e-100;
                     double x3_mirror = 0.;
                     double h1_mirror = (ICETHICK - FIRNDEPTH) + posnu[2];
                     double h2_mirror = FIRNDEPTH;
                     double nconst_mirror = n2 * n2 / n1 / n1; //here n2=NFIRN, n1=NICE
                     double x2_mirror = nconst_mirror;
                     double deltax_mirror;
//                     if (ST_TYPE == 4)
                        deltax_mirror = sqrt(Square(posnu[0] - MirrorATCoordinates8[WhichMirrorAntenna][0]) + Square(posnu[1] - MirrorATCoordinates8[WhichMirrorAntenna][1]));

//                     else {
//                        deltax_mirror = 0;
//                        cout << "Wrong ST_TYPE" << endl;
//                     }

                     if (deltax_mirror == 0) {
                        theta1_mirror = 0.;
                        theta2_mirror = 0.;
                        hy1_mirror = h1_mirror;
                        hy2_mirror = h2_mirror;
                     } else {
                        do {
                           x3_mirror = (x1_mirror + x2_mirror) / 2;
                           if (fx(x1_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)*fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror) < 0)
                              x2_mirror = x3_mirror;
                           else x1_mirror = x3_mirror;

                        } while (fabs(fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)) > 0.0001);

                        theta1_mirror = asin(sqrt(x3_mirror)); //incident angle in the ice,less than 90 degrees
                        theta2_mirror = asin(n1 * sin(theta1_mirror) / n2); //refraction angle in the firn, less than 90 degrees
                        hy1_mirror = h1_mirror / cos(theta1_mirror);
                        hy2_mirror = h2_mirror / cos(theta2_mirror);
                     }

                     theta_nposnu2MirrorAT = PI - theta1_mirror; //interacting in the ice, the theta angle close to the interaction point
                     d_posnu2MirrorAT = hy1_mirror + hy2_mirror;

                     double original_nsignal_mirror[3];
                     original_nsignal_mirror[0] = sin(theta_nposnu2MirrorAT) * cos(phi_nposnu2MirrorAT);
                     original_nsignal_mirror[1] = sin(theta_nposnu2MirrorAT) * sin(phi_nposnu2MirrorAT);
                     original_nsignal_mirror[2] = cos(theta_nposnu2MirrorAT);
                     //cout<<"whether the signal is downward"<<original_nsignal_mirror[2]<<endl;

                     viewangle_mirror = Angle(nnu, original_nsignal_mirror);

//start of new insert
//KD introduced here on 09/01 and 07/26/11 to get original modified polarization
                     GetPolarization(nnu, original_nsignal_mirror, n_pol, n_Bfield);


                     b1.original_nsignal_mirror[0] = original_nsignal_mirror[0];
                     b1.original_nsignal_mirror[1] = original_nsignal_mirror[1];
                     b1.original_nsignal_mirror[2] = original_nsignal_mirror[2];
                     b1.Polarization_R_original[0] = n_pol[0];
                     b1.Polarization_R_original[1] = n_pol[1];
                     b1.Polarization_R_original[2] = n_pol[2];
                     SphAngles(n_pol, sphangles);
                     b1.polariz_phi_R_original = sphangles[0] * RAD2DEG;
                     b1.polariz_theta_R_original = sphangles[1] * RAD2DEG;

                     /*
                     cout<<" >in ice "<<inu<<" "<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P1:"
                        <<original_nsignal_mirror[0]<<","<<original_nsignal_mirror[1]<<","<<original_nsignal_mirror[2]<<" angles:"
                        <<theta_nposnu2MirrorAT*RAD2DEG<<","<<phi_nposnu2MirrorAT*RAD2DEG<<">"<<viewangle_mirror*RAD2DEG<<" E1:"
                        <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<"="<<sphangles[1]*RAD2DEG<<","<<sphangles[0]*RAD2DEG
                        <<endl;  //theta2_mirror and  theta_nposnu2MirrorAT are supplementary


                        //after Snell's? //I think these 3 lines are redundant. maybe I just included them for cout that follows
                           nsignal_mirror_atAT[0]=sin(theta2_mirror)*cos(phi_nposnu2MirrorAT);
                           nsignal_mirror_atAT[1]=sin(theta2_mirror)*sin(phi_nposnu2MirrorAT);
                           nsignal_mirror_atAT[2]=-cos(theta2_mirror);

                     cout<<" >after firn "<<inu<<" nnu:"<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P2:"
                        <<nsignal_mirror_atAT[0]<<","<<nsignal_mirror_atAT[1]<<","<<nsignal_mirror_atAT[2]<<" angles:"
                        <<theta2_mirror*RAD2DEG<<","<<phi_nposnu2MirrorAT*RAD2DEG<<" E1 tree:"
                        <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]
                        <<endl;
                     */


//cout<<"inu "<<inu<<" R incident "<<(PI-theta_nposnu2MirrorAT)*RAD2DEG<<" transmitted "<<theta2_mirror*RAD2DEG<<endl;


                     if (DEPTH_DEPENDENT_N) {
                        GetRotation(PI - theta_nposnu2MirrorAT, theta2_mirror, original_nsignal_mirror, n_pol, n_pol_out);

                        /*
                        cout<<" >after rot "<<inu<<" nnu:"<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P1:"
                           <<original_nsignal_mirror[0]<<","<<original_nsignal_mirror[1]<<","<<original_nsignal_mirror[2]<<" Snells angles:"
                           <<theta2_mirror*RAD2DEG<<"<->"<<(PI-theta_nposnu2MirrorAT)*RAD2DEG<<" E1:"
                           <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<" E2:"
                           <<n_pol_out[0]<<","<<n_pol_out[1]<<","<<n_pol_out[2]
                           <<endl;
                        */

                     } else {
//GetRotation(PI-theta_nposnu2MirrorAT, theta2_mirror, original_nsignal_mirror, n_pol, n_pol_out);
// use above if we want to use rotation for non-graded binary indices
                        GetFresnel(PI - theta_nposnu2MirrorAT, theta2_mirror, original_nsignal_mirror, nsignal_mirror_atAT, n_pol, n_pol_out); //NOTE: we modify the incident angle definition
                     }

                     SphAngles(n_pol_out, sphangles);
                     b1.Fresnel_Pol_R[0] = n_pol_out[0];
                     b1.Fresnel_Pol_R[1] = n_pol_out[1];
                     b1.Fresnel_Pol_R[2] = n_pol_out[2];
                     b1.Fresnel_Pol_R_phi = sphangles[0] * RAD2DEG;
                     b1.Fresnel_Pol_R_theta = sphangles[1] * RAD2DEG;

                     /*
                     cout<<"Fresnel tree: "<<n_pol_out[0]<<","<<n_pol_out[1]<<","<<n_pol_out[2]<<"="
                              <<sphangles[1] * RAD2DEG<<","<<sphangles[0] * RAD2DEG
                              <<endl;
                     */

//equating these for populating trees below
                     n_pol[0] = n_pol_out[0];
                     n_pol[1] = n_pol_out[1];
                     n_pol[2] = n_pol_out[2];
//end of new insert
                  }//end of if(<400m), interact in the ice


                  else { //reflected events happen in the firn
                     double x1_mirror = 1.e-100;
                     double x3_mirror = 0.;
                     double h1_mirror = 2 * (ICETHICK - FIRNDEPTH);
                     double h2_mirror = FIRNDEPTH + (posnu[2] - (ICETHICK - FIRNDEPTH));
//double nconst_mirror=n2*n2/n1/n1;//this is wrong because in this case n1=n2=NFIRN since the interaction point is in the firn
//KD 09/06/11: NO FW, I think we should take into account the fact that interaction at at a depth and index changes
                     double nconst_mirror = NFIRN * NFIRN / NICE / NICE; //for FIRN case
                     double x2_mirror = nconst_mirror;
                     double deltax_mirror;
//                     if (ST_TYPE == 4)
                        deltax_mirror = sqrt(Square(posnu[0] - MirrorATCoordinates8[WhichMirrorAntenna][0]) + Square(posnu[1] - MirrorATCoordinates8[WhichMirrorAntenna][1]));
//                     else {
//                        deltax_mirror = 0;
//                        cout << "Wrong ST_TYPE" << endl;
//                     }

                     if (deltax_mirror == 0) {
                        theta1_mirror = 0.;
                        theta2_mirror = 0.;
                        hy1_mirror = h1_mirror;
                        hy2_mirror = h2_mirror;
                     } else {
                        do {
                           x3_mirror = (x1_mirror + x2_mirror) / 2;
                           if (fx(x1_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)*fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror) < 0)
                              x2_mirror = x3_mirror;
                           else x1_mirror = x3_mirror;

                        } while (fabs(fx(x3_mirror, h1_mirror, h2_mirror, nconst_mirror, deltax_mirror)) > 0.00001);

                        theta1_mirror = asin(sqrt(x3_mirror)); //incident angle in the ice,less than 90 degrees
//KD: 09/6/11 need to rethink these exact calculations
//theta2_mirror=asin(n1*sin(theta1_mirror)/n2);//refraction angle in the firn, less than 90 degrees
                        theta2_mirror = asin(NICE * sin(theta1_mirror) / NFIRN);
                        hy1_mirror = h1_mirror / cos(theta1_mirror);
                        hy2_mirror = h2_mirror / cos(theta2_mirror);

                        // cout<<hy1_mirror<<" "<<hy2_mirror<<endl;

                     }
                     theta_nposnu2MirrorAT = PI - theta2_mirror; //the theta angle at the mirror antennas, also the theta angle of the signal at the interaction point since interating in the firn
                     d_posnu2MirrorAT = hy1_mirror + hy2_mirror;

                     double original_nsignal_mirror[3];
                     original_nsignal_mirror[0] = sin(theta_nposnu2MirrorAT) * cos(phi_nposnu2MirrorAT);
                     original_nsignal_mirror[1] = sin(theta_nposnu2MirrorAT) * sin(phi_nposnu2MirrorAT);
                     original_nsignal_mirror[2] = cos(theta_nposnu2MirrorAT);
                     //cout<<"whether the signal is downward"<<original_nsignal_mirror[2]<<endl;

                     viewangle_mirror = Angle(nnu, original_nsignal_mirror);
//KD introduced here to get original modified polarization
                     GetPolarization(nnu, original_nsignal_mirror, n_pol, n_Bfield);
//previously it was getting polarization for LOS ray by default; this is now slightly better
                     b1.original_nsignal_mirror[0] = original_nsignal_mirror[0];
                     b1.original_nsignal_mirror[1] = original_nsignal_mirror[1];
                     b1.original_nsignal_mirror[2] = original_nsignal_mirror[2];
                     b1.Polarization_R_original[0] = n_pol[0];
                     b1.Polarization_R_original[1] = n_pol[1];
                     b1.Polarization_R_original[2] = n_pol[2];
                     SphAngles(n_pol, sphangles);
                     b1.polariz_phi_R_original = sphangles[0] * RAD2DEG;
                     b1.polariz_theta_R_original = sphangles[1] * RAD2DEG;


                     /*
                     cout<<" >at vertex "<<inu<<" "<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P1:"
                        <<original_nsignal_mirror[0]<<","<<original_nsignal_mirror[1]<<","<<original_nsignal_mirror[2]<<" P1 angles:"
                        <<theta_nposnu2MirrorAT*RAD2DEG<<","<<phi_nposnu2MirrorAT*RAD2DEG<<">"<<viewangle_mirror*RAD2DEG<<" E1:"
                        <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<" = " <<sphangles[1]*RAD2DEG<<","<<sphangles[0]*RAD2DEG<<endl
                        <<"theta1_mirror :"<<theta1_mirror*RAD2DEG<<"indices :"<<n1<<" "<<n2<<" N:"<<NICE<<" "<<NFIRN
                        <<endl;  //theta2_mirror and  theta_nposnu2MirrorAT are supplementary
                     //after Snell's? //I think these 3 lines are redundant. maybe I just included them for cout that follows
                           nsignal_mirror_atAT[0]=sin(theta2_mirror)*cos(phi_nposnu2MirrorAT);
                           nsignal_mirror_atAT[1]=sin(theta2_mirror)*sin(phi_nposnu2MirrorAT);
                           nsignal_mirror_atAT[2]=-cos(theta2_mirror);

                     cout<<" >after firn "<<inu<<" nnu:"<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P2:"
                        <<nsignal_mirror_atAT[0]<<","<<nsignal_mirror_atAT[1]<<","<<nsignal_mirror_atAT[2]<<" P2  angles:"
                        <<theta2_mirror*RAD2DEG<<","<<phi_nposnu2MirrorAT*RAD2DEG<<" E1 tree:"
                        <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]
                        <<endl;

                     */

//here we introduce a similar process to the ice vertex as above
                     if (DEPTH_DEPENDENT_N) {
                        GetRotation(PI - theta_nposnu2MirrorAT, theta2_mirror, original_nsignal_mirror, n_pol, n_pol_out);
                        /*
                        cout<<" >after rot "<<inu<<" nnu:"<<nnu[0]<<","<<nnu[1]<<","<<nnu[2]<<" ; P1:"
                           <<original_nsignal_mirror[0]<<","<<original_nsignal_mirror[1]<<","<<original_nsignal_mirror[2]<<" Snells angles:"
                           <<theta2_mirror*RAD2DEG<<"<->"<<(PI-theta_nposnu2MirrorAT)*RAD2DEG<<" E1:"
                           <<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<" E2:"
                           <<n_pol_out[0]<<","<<n_pol_out[1]<<","<<n_pol_out[2]
                           <<endl;
                        */

                     } else {
//GetRotation(PI-theta_nposnu2MirrorAT, theta2_mirror, original_nsignal_mirror, n_pol, n_pol_out);
// use above if we want to use rotation for non-graded binary indices
                        GetFresnel(PI - theta_nposnu2MirrorAT, theta2_mirror, original_nsignal_mirror, nsignal_mirror_atAT, n_pol, n_pol_out); //NOTE: we modify the incident angle definition
                     }

                     SphAngles(n_pol_out, sphangles);
                     b1.Fresnel_Pol_R[0] = n_pol_out[0];
                     b1.Fresnel_Pol_R[1] = n_pol_out[1];
                     b1.Fresnel_Pol_R[2] = n_pol_out[2];
                     b1.Fresnel_Pol_R_phi = sphangles[0] * RAD2DEG;
                     b1.Fresnel_Pol_R_theta = sphangles[1] * RAD2DEG;

                     /*
                     cout<<"Fresnel tree: "<<n_pol_out[0]<<","<<n_pol_out[1]<<","<<n_pol_out[2]<<"="
                              <<sphangles[1] * RAD2DEG<<","<<sphangles[0] * RAD2DEG
                              <<endl;

                     */
//equating these for populating trees below
                     n_pol[0] = n_pol_out[0];
                     n_pol[1] = n_pol_out[1];
                     n_pol[2] = n_pol_out[2];
                  }//end of reflected events in firn (note order is different in direct)

                  nsignal_mirror_atAT[0] = sin(theta2_mirror) * cos(phi_nposnu2MirrorAT);
                  nsignal_mirror_atAT[1] = sin(theta2_mirror) * sin(phi_nposnu2MirrorAT);
                  nsignal_mirror_atAT[2] = -cos(theta2_mirror);

                  //08/30/11 As specified above, this following procedure to obtain polarization in firn in NOT correct
                  // GetPolarization(nnu,nsignal_mirror_atAT,n_pol,n_Bfield);

               }//END of if(FIRN)





               if (PLANEWAVE) {
//HAVE TO REVIEW THIS CALCULATION FOR REFLECTED

//calculating distance and time from AT to STcenter, using distance of pt from a plane, where plane is defined by its normal
// =|ax2 + by2 + cy2 +d| / sqrt(a^2+b^2+c^2)
                  double d_AT2ST_mirror_denominator = (nsignal_atMirrorST[0] * MirrorATCoordinates8[WhichMirrorAntenna][0] + nsignal_atMirrorST[1] * MirrorATCoordinates8[WhichMirrorAntenna][1] +
                                                       nsignal_atMirrorST[2] * MirrorATCoordinates8[WhichMirrorAntenna][2] + (-nsignal_atMirrorST[0] * MirrorATCoordinate[0] -
                                                             -nsignal_atMirrorST[1] * MirrorATCoordinate[1] - nsignal_atMirrorST[2] * MirrorATCoordinate[2]));
                  double d_AT2ST_mirror_numerator = sqrt(Square(nsignal_atMirrorST[0]) + Square(nsignal_atMirrorST[1]) + Square(nsignal_atMirrorST[2]));
                  double d_AT2ST_mirror = d_AT2ST_mirror_denominator / d_AT2ST_mirror_numerator ; //why do I have it inverted??

                  double time_AT2ST_mirror = 0; //time offset compared to center of ST
                  time_AT2ST_mirror = d_AT2ST_mirror * NICE / C;
                  if (FIRN)
                     time_AT2ST_mirror = d_AT2ST_mirror * NFIRN / C; //station is always in firn, right?

                  b1.my_time_AT2ST_mirror[WhichMirrorAntenna] = time_AT2ST_mirror;
                  b1.nsignal_mirror_i_atAT[WhichMirrorAntenna] =  nsignal_mirror_atAT[0];
                  b1.nsignal_mirror_j_atAT[WhichMirrorAntenna] =  nsignal_mirror_atAT[1];
                  b1.nsignal_mirror_k_atAT[WhichMirrorAntenna] =  nsignal_mirror_atAT[2];
                  b1.nposnu2MirrorAT_i[WhichMirrorAntenna] =  nposnu2MirrorAT[0];
                  b1.nposnu2MirrorAT_j[WhichMirrorAntenna] =  nposnu2MirrorAT[1];
                  b1.nposnu2MirrorAT_k[WhichMirrorAntenna] =  nposnu2MirrorAT[2];
               }//END of if(PLANEWAVE)


               b1.Polarization_R[0] = n_pol[0];
               b1.Polarization_R[1] = n_pol[1];
               b1.Polarization_R[2] = n_pol[2];

               SphAngles(n_pol, sphangles);
               b1.polariz_phi_R = sphangles[0] * RAD2DEG;
               b1.polariz_theta_R = sphangles[1] * RAD2DEG;

               /*
               cout<<"E2 tree:"<<n_pol[0]<<","<<n_pol[1]<<","<<n_pol[2]<<"="
                     <<sphangles[1] * RAD2DEG<<","<<sphangles[0] * RAD2DEG
                     <<endl;
               */


               Zero(vmmhz, NFREQ);
               Zero(vmmhz1m, NFREQ);
               Zero(vmmhz1m_unattened, NFREQ);
               double sum_vmmhz_mirror = 0.;
               double sum_vmmhz1m_unattened_mirror = 0;
               double sum_vmmhz1m_unattened_mirror_max = 0;
               double volt_mirror_max = 0.;

               double   Esumvmmhz1m_unattened_mirror = 0.;

               for (int i = 0; i < NFREQ; i++) {
                  vmmhz1m_unattened[i] = GetVmMHz1m(pnu, freq[i], X0ICE, ECICE, n1, AEX_ICE);
//-------------------------------
                  b7.Evmmhz1m_unattened_mirror[i] = vmmhz1m_unattened[i];
//-------------------------------
                  vmmhz1m_unattened_max[i] = vmmhz1m_unattened[i];
                  Esumvmmhz1m_unattened_mirror += FREQ_BIN * vmmhz1m_unattened[i]; //same as for direct as this variable, WELL MAYBE NOT, if it escapes the direct loop.
               }
//-------------------------------
               b7.Esumvmmhz1m_unattened_mirror = Esumvmmhz1m_unattened_mirror;
//-------------------------------

               if (SECKEL == 0) {
                  vmmhz_max = VmMHz(vmmhz1m_max, d_posnu2MirrorAT);
                  if (ATTEN_FREQ) {
                     for (int i = 0; i < NFREQ; i++) {
                        if (freq[i] <= 122) {
                           attenlength_down_freq = attenlength_down *= (700. - 250.*(freq[i] - 50.) / 72.) / 450.; //put 300 here artificially

                           vmmhz_max_freq[i] = VmMHz_attenuated(d_posnu2MirrorAT, vmmhz_max, attenlength_down_freq);
                        } else
                           vmmhz_max_freq[i] = VmMHz_attenuated(d_posnu2MirrorAT, vmmhz_max, attenlength_down);

                     }


                  }

                  vmmhz_max = VmMHz_attenuated(d_posnu2MirrorAT, vmmhz_max, attenlength_down);
                  if (ATTEN_FREQ) {
//                     if (ST_TYPE == 4) {
                        if (vmmhz_max_freq[0]*heff_max_LPA * BW < NV) {
                           continue;
                        }
//                     }

                     GetVmMHz_freq(vmmhz_max_freq, vmmhz1m_max, pnu, vmmhz);

                  } else {
//                     if (ST_TYPE == 4) {
                        if (vmmhz_max * heff_max_LPA * BW < NV) {
                           continue;
                        }
//                     }
                     GetVmMHz(vmmhz_max, vmmhz1m_max, pnu, vmmhz);
//KD these are different now because mirror distances and attenuations come in
                     double Esum_vmmhz_mirror = 0.;
                     for (int i = 0; i < NFREQ; i++) {
//-------------------------------
                        b7.Evmmhz_mirror[WhichMirrorAntenna][i] = vmmhz[i];
//-------------------------------
                        Esum_vmmhz_mirror += FREQ_BIN * vmmhz[i];
                     }

//-------------------------------
                     b7.Esum_vmmhz_mirror[WhichMirrorAntenna] = Esum_vmmhz_mirror;//=(attenfactor*Esumvmmhz1m_unattened/dis)  //-------------------------------

                  }
               }

               for (int i = 0; i < NFREQ; i++) {
                  if (SECKEL == 0) {

                     vmmhz[i] *= sqrt(REFLECT_RATE); //50% power reflected,3dB
                     // cout<<vmmhz[i]<<endl;
                     // cout<<endl;
                     //  vmmhz[i]*=0.32;//10% power reflected, 10dB
                     // vmmhz[i]*=0.1; //1%power reflected,20dB
                     TaperVmMHz(viewangle_mirror, deltheta_em[i], deltheta_had[i], emfrac, hadfrac, vmmhz[i], vmmhz_em[i], vmmhz_had[i]);
//-------------------------------
                     b7.Evmmhz_taper_mirror[WhichMirrorAntenna][i] = vmmhz[i];
                     b7.Evmmhz_em_mirror[WhichMirrorAntenna][i] = vmmhz_em[i];
                     b7.Evmmhz_had_mirror[WhichMirrorAntenna][i] = vmmhz_had[i];
//-------------------------------
                     TaperVmMHz(viewangle_mirror, deltheta_em[i], deltheta_had[i], emfrac, hadfrac, vmmhz1m_unattened[i], vmmhz_em[i], vmmhz_had[i]);
//b1.Efield_atposnu_mirror_my[WhichMirrorAntenna][i] = vmmhz1m_unattened[i];
//-------------------------------
                     b7.Efield_atposnu_mirror_my[WhichMirrorAntenna][i] = vmmhz1m_unattened[i];
//-------------------------------
                     TaperVmMHz(changle, deltheta_em[i], deltheta_had[i], emfrac, hadfrac, vmmhz1m_unattened_max[i], vmmhz_em[i], vmmhz_had[i]);
                  }

                  if (SECKEL == 1) {
                     vmmhz1m[i] = VmMHz1m(viewangle_mirror, freq[i], emfrac, hadfrac, em_shower_length, had_shower_length);
                     vmmhz[i] = VmMHz(vmmhz1m[i], d_posnu2MirrorAT);
                     vmmhz[i] = 0.7 * VmMHz_attenuated(d_posnu2MirrorAT, vmmhz[i], attenlength_down);
                  }


                  sum_vmmhz_mirror += (FREQ_BIN * vmmhz[i]);
                  sum_vmmhz1m_unattened_mirror += (FREQ_BIN * vmmhz1m_unattened[i]);
                  sum_vmmhz1m_unattened_mirror_max += (FREQ_BIN * vmmhz1m_unattened_max[i]);

               }//end of freq for loop
//-------------------------------
               b7.Esum_vmmhz_taper_mirror[WhichMirrorAntenna] = sum_vmmhz_mirror;
//-------------------------------

               abs_time_mirror = (hy1_mirror * NICE + hy2_mirror * NFIRN) / C;

//cout<<"|"<<inu<<"(REF):"<<iRow_oncone_mirror.at(WhichMirrorStation)<<","<<iCol_oncone_mirror.at(WhichMirrorStation)<<","<<NV;

//               if (ST_TYPE == 4) {
                  volt_mirror_max = sum_vmmhz_mirror * heff_max_LPA;
                  if (volt_mirror_max < NV) {
                     b1.iAT_mirror[WhichMirrorAntenna]      = 0;
                     //b1.kviewangle_mirror[WhichMirrorAntenna]= 0;
                     //b1.volt_mirror_max[WhichMirrorAntenna]  = 0;
                     b1.volt_LPA_mirror[WhichMirrorAntenna] = 0;
                     //b1.sum_vmmhz_mirror[WhichMirrorAntenna] = 0;
                     //b1.sum_vmmhz1m_unattened_mirror[WhichMirrorAntenna] = 0;
                     b1.hitangle_e_LPA_mirror[WhichMirrorAntenna] = 0;
                     b1.hitangle_h_LPA_mirror[WhichMirrorAntenna] = 0;
                     b1.totaltime_mirror[WhichMirrorAntenna]      = 0;

                     continue;
                  }
//               }

               //int count_channels_mirror=0; // KD: a relic of ST_TYPE 3?


//Replace with new antenna model               if (ST_TYPE == 4) {

/*  //OldShelfMC
                  if (FIRN)
                     GetMirrorHitAngle_LPA(WhichMirrorAntenna, N_Ant_perST, nsignal_mirror_atAT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA);
                  else
                     GetMirrorHitAngle_LPA(WhichMirrorAntenna, N_Ant_perST, nposnu2MirrorAT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA);
                  double volt_LPA_mirror = 0; //volts of log periodic antenna
                  double volt_LPA_mirror_preNoise = 0;
                  term_LPA = 0; //zero term_LPA
                  term_LPA_e = 0;
                  term_LPA_h = 0;

//KD are these different for direct and reflected? Yep
                  b1.e_component_LPA_mirror[WhichMirrorAntenna] = e_component_LPA;
                  b1.h_component_LPA_mirror[WhichMirrorAntenna] = h_component_LPA;
*/


	       //antenna orientation vectors to pass into Heff function
	          double n_boresight_mirror[3];
		  double n_eplane_mirror[3];
		  for (int i =0; i<3; i++){
		    n_boresight_mirror[i]=MirrorAnt_n_boresight[WhichMirrorAntenna][i];
		    n_eplane_mirror[i]=MirrorAnt_n_eplane[WhichMirrorAntenna][i];
		  }

		  //Give an e and h plane component for output tree
		  if (FIRN)
		    GetHitAngle(n_boresight_mirror, n_eplane_mirror, nsignal_mirror_atAT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA);
                  else
		    GetHitAngle(n_boresight_mirror, n_eplane_mirror, nposnu2MirrorAT, n_pol, hitangle_e_LPA, hitangle_h_LPA, e_component_LPA, h_component_LPA);
		  
                  double volt_LPA_mirror = 0; //volts of log periodic antenna
                  double volt_LPA_mirror_preNoise = 0;
                  term_LPA = 0; //zero term_LPA
                  term_LPA_e = 0;
                  term_LPA_h = 0;

//KD are these different for direct and reflected? Yep
                  b1.e_component_LPA_mirror[WhichMirrorAntenna] = e_component_LPA;
                  b1.h_component_LPA_mirror[WhichMirrorAntenna] = h_component_LPA;


                  for (int i = 0; i < NFREQ; i++) {

//term_LPA=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(gainv,freq[i]*1.E6)*
//sqrt(       pow(e_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[0][i])*(hitangle_e_LPA/flare[0][i]))*exp(-2*ALOG2*(hitangle_h_LPA/flare[1][i])*(hitangle_h_LPA/flare[1][i]))   ,2)     );

/* //2-23-2017 moved to GetHeff()
                     term_LPA = vmmhz[i] * FREQ_BIN * 0.5 * GaintoHeight(gainv, freq[i] * 1.E6) *
                                sqrt((pow(e_component_LPA * exp(-2 * ALOG2 * (hitangle_e_LPA / flare[0][i]) * (hitangle_e_LPA / flare[0][i])), 2)  +   pow(e_component_LPA * exp(-2 * ALOG2 * (hitangle_h_LPA / flare[1][i]) * (hitangle_h_LPA / flare[1][i])), 2)) / 2);
*/

		    if (FIRN)
		      term_LPA = vmmhz[i] * FREQ_BIN * 0.5 * GetHeff(AntType[WhichMirrorAntenna], freq[i],n_boresight_mirror,n_eplane_mirror, nsignal_mirror_atAT, n_pol);
		    else
		      term_LPA = vmmhz[i] * FREQ_BIN * 0.5 * GetHeff(AntType[WhichMirrorAntenna], freq[i],n_boresight_mirror,n_eplane_mirror, nposnu2MirrorAT, n_pol);

                     /*//works for gain tests
                     term_LPA=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(gainv,freq[i]*1.E6)*sqrt(pow(e_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[0][i])*(hitangle_e_LPA/flare[0][i])),2)     +
                     0.01*pow(h_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[1][i])*(hitangle_e_LPA/flare[1][i])),2));
                     */


                     /*
                     //term_LPA_h=term_LPA_e;
                     term_LPA_h=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(gainv,freq[i]*1.E6)*sqrt(pow(h_component_LPA*exp(-2*ALOG2*(hitangle_h_LPA/flare[1][i])*(hitangle_h_LPA/flare[1][i])),2) + 0.01*pow(e_component_LPA*exp(-2*ALOG2*(hitangle_h_LPA/flare[0][i])*(hitangle_h_LPA/flare[0][i])),2));
                     term_LPA = 0.5*(term_LPA_e + term_LPA_h);                   */

                     /*

                     // factor of 2 should be in exponent

                     //&&&&&&&& OLD gain equation
                     term_LPA=vmmhz[i]*FREQ_BIN*0.5*GaintoHeight(GetGainV(freq[i]*1.E6),freq[i]*1.E6)*sqrt(pow(e_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[0][i])*(hitangle_e_LPA/flare[0][i])),2)
                                                 + 0.01*pow(h_component_LPA*exp(-2*ALOG2*(hitangle_e_LPA/flare[1][i])*(hitangle_e_LPA/flare[1][i])),2));;
                     */
//KD heff should be the same for direct
//b1.heff_my[WhichAntenna][i] = GaintoHeight(gainv,freq[i]*1.E6);
//-------------------------------
//b7.heff_my[i]   = GaintoHeight(gainv,freq[i]*1.E6);
//-------------------------------

                     //-----------------------
                     b1.vmmhz_mirror_my[WhichMirrorAntenna][i] = vmmhz[i];
                     b1.term_LPA_mirror_my[WhichMirrorAntenna][i] = term_LPA;
//-------------------------------
                     b7.vmmhz_mirror_my[WhichMirrorAntenna][i] = vmmhz[i];
                     b7.term_LPA_mirror_my[WhichMirrorAntenna][i] = term_LPA;
//-------------------------------
                     volt_LPA_mirror += term_LPA;

                  }

//if (eshower_em > 1.E10){
//cout<<"inu= "<<inu<<" pnu, em-, had-shower= ("<<pnu<<","<<eshower_em<<","<<eshower_had<<") volt -max_M, -LPA_M= ("<<volt_mirror_max<<","<<volt_LPA_mirror<<") changle, viewangle_M= ("<<changle*RAD2DEG<<","<<viewangle_mirror*RAD2DEG<<") sum_vmmhz _M,unattend_M, unattened_max_M = ("<<sum_vmmhz_mirror<<","<<sum_vmmhz1m_unattened_mirror<<","<<sum_vmmhz1m_unattened_mirror_max<<")"<<endl;
//}


                  if (SIGNAL_FLUCT) {  //KD 10/11/2010 adding time precision fluctuation here
                     volt_LPA_mirror_preNoise = volt_LPA_mirror; //gets populated if Noise present
                     volt_LPA_mirror += Rand3.Gaus(0., VNOISE);
                     abs_time_mirror += Rand3.Gaus(0., TIMEPRECISION);
                  }


//cout<<"||"<<inu<<"(REF):"<<iRow_oncone_mirror.at(WhichMirrorStation)<<","<<iCol_oncone_mirror.at(WhichMirrorStation)<<","<<NV<<endl;

                  if (volt_LPA_mirror >= NV) {
                     sum_triggeredAT_mirror++;
                     iAT_mirror_L1.push_back(WhichMirrorAntenna);
                     //cout<<WhichMirrorAntenna<<endl;
                     //cout<<hy1_mirror<<" "<<hy2_mirror<<endl;
                     path_inice_mirror_L1.push_back(hy1_mirror);
                     path_infirn_mirror_L1.push_back(hy2_mirror);
                     totaltime_mirror_L1.push_back(abs_time_mirror);
                     Efield_atposnu_mirror_L1.push_back(sum_vmmhz1m_unattened_mirror);
                     volte_allAT_mirror_L1.push_back(volt_LPA_mirror);
                     volth_allAT_mirror_L1.push_back(0);
                     AT_posi_x_mirror_L1.push_back(MirrorATCoordinates8[WhichMirrorAntenna][0]);
                     AT_posi_y_mirror_L1.push_back(MirrorATCoordinates8[WhichMirrorAntenna][1]);
                     AT_posi_z_mirror_L1.push_back(MirrorATCoordinates8[WhichMirrorAntenna][2]);
                     Dir_shower_x_mirror_L1.push_back(nsignal_mirror_atAT[0]);
                     Dir_shower_y_mirror_L1.push_back(nsignal_mirror_atAT[1]);
                     Dir_shower_z_mirror_L1.push_back(nsignal_mirror_atAT[2]);
                     theta_shower_atposnu_mirror_L1.push_back(theta1_mirror);
                     hitangle_e_LPA_mirror_L1.push_back(hitangle_e_LPA * RAD2DEG);
                     hitangle_h_LPA_mirror_L1.push_back(hitangle_h_LPA * RAD2DEG);
                  }

//b1.kviewangle_mirror[WhichMirrorAntenna]= viewangle_mirror*RAD2DEG;
//b1.volt_mirror_max[WhichMirrorAntenna]  = volt_mirror_max;
                  b1.volt_LPA_mirror[WhichMirrorAntenna] = volt_LPA_mirror;
                  b1.volt_LPA_mirror_preNoise[WhichMirrorAntenna] = volt_LPA_mirror_preNoise;
//b1.sum_vmmhz_mirror[WhichMirrorAntenna] = sum_vmmhz_mirror;
//b1.sum_vmmhz1m_unattened_mirror[WhichMirrorAntenna] = sum_vmmhz1m_unattened_mirror;
//b1.hitangle_e_LPA_mirror[WhichMirrorAntenna]     = hitangle_e_LPA;

//replace with new antenna model               }

               b1.sum_vmmhz1m_unattened_mirror_max       = sum_vmmhz1m_unattened_mirror_max;
               b1.hitangle_e_LPA_mirror[WhichMirrorAntenna] = hitangle_e_LPA * RAD2DEG;
               b1.hitangle_h_LPA_mirror[WhichMirrorAntenna] = hitangle_h_LPA * RAD2DEG;
               b1.theta_my_signalAT_mirror[WhichMirrorAntenna] = theta2_mirror * RAD2DEG;
               b1.phi_my_signalAT_mirror[WhichMirrorAntenna]   = phi_nposnu2MirrorAT * RAD2DEG;
               b1.totaltime_mirror[WhichMirrorAntenna]      = abs_time_mirror;
//b1.volt_LPA_mirror[WhichMirrorAntenna]  = volte_allAT_mirror_L1.at(volte_allAT_mirror_L1.size()-1);
            }//end of  antennas FOR loop per station

            //cout<<"4.4th place"<<double(clock()-init_time)/CLOCKS_PER_SEC<<endl;

//++cout <<"KD  : Having populated L1_mirror arrays, I can now output L1_mirror details:"<<endl;

//++for (int i=0; i<iAT_mirror_L1.size(); i++) //KD: looping over triggered antenna
//++{
//++cout<<"KD21a: "<<"iAT_mirror_L1[i]:"<< iAT_mirror_L1[i] << "  AT mirror position L1[i]=("<< AT_posi_x_mirror_L1[i] << ","<< AT_posi_y_mirror_L1[i] << "," << AT_posi_z_mirror_L1[i] <<")  volt_LPA_mirror="<< volte_allAT_mirror_L1[i]<< " totaltime_mirror_L1="<<totaltime_mirror_L1[i]<<endl;
//"  Efield_atposnu_mirror_L1[i]="<< Efield_atposnu_mirror_L1[i] << endl;
//++}

//+cout<<"largest time difference for L1_mirror = maxi - mini = "<<maxi<<"-"<<mini<<" = "<<maxi-mini<<endl;

//++cout<<"KD22b: sum_triggeredAT_mirror: " << sum_triggeredAT_mirror << endl;
//cout<<"------------------------------------------------------------------------------"<<endl;


//++cout<<"KD  : If that sum_t'AT_mirror is <2, continue; else more print statements..."<<endl;
            if (sum_triggeredAT_mirror < N_Ant_Trigger - 1) //KD: new variable
//       if(sum_triggeredAT_mirror<2)//if the mirror station is not triggered
               continue;
            else {
//               if (ST_TYPE == 4) {
                  if (sum_triggeredAT_mirror == N_Ant_Trigger - 1) //KD : new variable
//          if(sum_triggeredAT_mirror==2)//at least 3 out of 8
                     continue;
//               }
//++cout<<"KD  : ...where we continue because sum_t'_AT_mirror isn't 2 but is: "<<sum_triggeredAT_mirror<<endl;

               if (0) { // activate for studies of different scenario
//-------------------------------------------------------------
//-----STUART's question---------------------------------------
                  int even_Mant = 0; //resetting variable
                  int odd_Mant = 0; //resetting variable

                  for (uint WhichMAnt = 0; WhichMAnt < iAT_mirror_L1.size(); WhichMAnt++) {
//go through all triggered AT sequentially
                     //if(iAT_mirror_L1[WhichMAnt]==0||iAT_mirror_L1[WhichMAnt]==2
                     //||iAT_mirror_L1[WhichMAnt]==4||iAT_mirror_L1[WhichMAnt]==6)
                     if (iAT_mirror_L1[WhichMAnt] == 0) { //ignoring AT=2 AT=4
                        even_Mant++;//incrementing count for even antennas
                     }

                     if (iAT_mirror_L1[WhichMAnt] == 1 || iAT_mirror_L1[WhichMAnt] == 3 || iAT_mirror_L1[WhichMAnt] == 5) { //ignoring AT=5
                        odd_Mant++;//incrementing count for even antennas
                     }
                  }
//after going through all loops and getting even and odd triggered number of ant's
                  if (even_Mant < 1 && odd_Mant < 2) {
                     continue; //discard if not 1 AND 3, go out of loop and do not trigger station
                  }
//------------------------------------------------------------
               }

               sum_triggeredST_mirror++;
               N_TriggeredAT_perST_mirror_L2.push_back(sum_triggeredAT_mirror);
               Trig_Type_mirror_L2.push_back(-1);
               iRow_mirror_L2.push_back(iRow_oncone_mirror.at(WhichMirrorStation));
               iCol_mirror_L2.push_back(iCol_oncone_mirror.at(WhichMirrorStation));
               //  cout<<"line 3197"<<endl;
//++cout<<"KD  : The code has just updated 'sum_triggeredST_mirror' and filled first four L2_mirror vectors."<<endl;
               for (int i = 0; i < sum_triggeredAT_mirror; i++) {
                  // cout<<sum_triggeredAT_mirror<<"  "<<iAT_mirror_L1.size()<<endl;
                  //  cout<<"line 3201"<<endl;
                  //cout<<path_infirn_mirror_L1.at(i)<<endl;
                  iAT_mirror_L2.push_back(iAT_mirror_L1.at(i));
                  path_inice_mirror_L2.push_back(path_inice_mirror_L1.at(i));
                  path_infirn_mirror_L2.push_back(path_infirn_mirror_L1.at(i));
                  totaltime_mirror_L2.push_back(totaltime_mirror_L1.at(i));
                  //  cout<<"line 3206"<<endl;
                  Efield_atposnu_mirror_L2.push_back(Efield_atposnu_mirror_L1.at(i));
                  volte_allAT_mirror_L2.push_back(volte_allAT_mirror_L1.at(i));
                  volth_allAT_mirror_L2.push_back(volth_allAT_mirror_L1.at(i));
                  AT_posi_x_mirror_L2.push_back(AT_posi_x_mirror_L1.at(i));
                  AT_posi_y_mirror_L2.push_back(AT_posi_y_mirror_L1.at(i));
                  AT_posi_z_mirror_L2.push_back(AT_posi_z_mirror_L1.at(i));
                  //  cout<<"line 3213"<<endl;
                  Dir_shower_x_mirror_L2.push_back(Dir_shower_x_mirror_L1.at(i));
                  Dir_shower_y_mirror_L2.push_back(Dir_shower_y_mirror_L1.at(i));
                  Dir_shower_z_mirror_L2.push_back(Dir_shower_z_mirror_L1.at(i));
                  theta_shower_atposnu_mirror_L2.push_back(theta_shower_atposnu_mirror_L1.at(i));
//++cout<<"KD  : This repeats sum_t'AT_mirror times as I fill up remaining fourteen L2_mirror vectors.+"<<endl;
                  hitangle_e_LPA_mirror_L2.push_back(hitangle_e_LPA_mirror_L1.at(i));
                  hitangle_h_LPA_mirror_L2.push_back(hitangle_h_LPA_mirror_L1.at(i));
               }
               //  cout<<"line 3216"<<endl;

               double Efield_atposnu_mirror_L2_mymax  = 0;
               double volte_allAT_mirror_L2_mymax  = 0;
               double volte_allAT_mirror_L2_mymin  = 31.55;
               double hitangle_e_LPA_mirror_mymax  = 0;
               double hitangle_h_LPA_mirror_mymax  = 0;
               double totaltime_mirror_L2_mymax = 0;
               double totaltime_mirror_L2_mymin = 31.55; //arbirtrary max


               for (uint j = 0; j < iAT_mirror_L2.size(); j++) {
//cout<<"iAT_mirror = "<< iAT_mirror_L2.at(j)<< " Efield_atposnu_mirror= " << Efield_atposnu_mirror_L2.at(j) << " volt_e_mirror= "<< volte_allAT_mirror_L2.at(j)<<endl;
                  b1.iAT_mirror[j]        = iAT_mirror_L2[j] + 1;

                  if (Efield_atposnu_mirror_L2[j] > Efield_atposnu_mirror_L2_mymax)
                     Efield_atposnu_mirror_L2_mymax =  Efield_atposnu_mirror_L2[j];

                  if (volte_allAT_mirror_L2[j] >  volte_allAT_mirror_L2_mymax) {
                     volte_allAT_mirror_L2_mymax   = volte_allAT_mirror_L2[j];
                     hitangle_e_LPA_mirror_mymax   = hitangle_e_LPA_mirror_L2[j];
                     hitangle_h_LPA_mirror_mymax   = hitangle_h_LPA_mirror_L2[j];
                     iAT_mirror_index_mymax     = iAT_mirror_L2[j];
                  }

                  if (volte_allAT_mirror_L2[j] < volte_allAT_mirror_L2_mymin) {
                     volte_allAT_mirror_L2_mymin  = volte_allAT_mirror_L2[j];
                     iAT_mirror_index_mymin     = iAT_mirror_L2[j];
                  }

                  if (totaltime_mirror_L2[j] > totaltime_mirror_L2_mymax) {
                     totaltime_mirror_L2_mymax   = totaltime_mirror_L2[j];
                     iAT_mirror_index_timemax = iAT_mirror_L2[j];
                  }

                  if (totaltime_mirror_L2[j] < totaltime_mirror_L2_mymin) {
                     totaltime_mirror_L2_mymin = totaltime_mirror_L2[j];
                     iAT_mirror_index_timemin = iAT_mirror_L2[j];
                  }
               }

//cout<<"Efield_atposnu_mirror_L2_mymax= "<< Efield_atposnu_mirror_L2_mymax <<" volte_allAT_mirror_L2_mymax= "<<volte_allAT_mirror_L2_mymax<<endl;
               b1.Efield_atposnu_mirror_mymax      = Efield_atposnu_mirror_L2_mymax;
               b1.volte_allAT_mirror_mymax      = volte_allAT_mirror_L2_mymax;
               b1.iAT_mirror_index_mymax     = iAT_mirror_index_mymax;
               b1.iAT_mirror_index_mymin     = iAT_mirror_index_mymin;
               b1.hitangle_e_LPA_mirror_mymax      = hitangle_e_LPA_mirror_mymax;
               b1.hitangle_h_LPA_mirror_mymax      = hitangle_h_LPA_mirror_mymax;
               b1.totaltime_mirror_mymax     = totaltime_mirror_L2_mymax;
               b1.totaltime_mirror_mymin     = totaltime_mirror_L2_mymin;
               b1.iAT_mirror_index_timemax      = iAT_mirror_index_timemax;
               b1.iAT_mirror_index_timemin      = iAT_mirror_index_timemin;
            }

//++cout <<"KD  : I'm now outside the 'ELSE sum_AT_mirror<2' loop with sum_triggeredAT_mirror still: "<<sum_triggeredAT_mirror<<endl;
//++cout <<"KD  : I can now output L2_mirror details, that should be exact as L1_mirror above:"<<endl;
//++for (int i=0; i<sum_triggeredAT_mirror; i++) //KD: looping over triggered antenna
//++{
//++cout<<"KD22a: "<<"iAT_mirror_L2[i]:"<< iAT_mirror_L2[i] << "  AT mirror position L2[i]=("<< AT_posi_x_mirror_L2[i] << ","<< AT_posi_y_mirror_L2[i] << "," << AT_posi_z_mirror_L2[i] <<")  volt_LPA_mirror="<< volte_allAT_mirror_L2[i]<< " totaltime_mirror_L2="<<totaltime_mirror_L2[i]<<endl;
//"  Efield_atposnu_L2[i]="<< Efield_atposnu_mirror_L2[i] << endl;
//++}
//++cout<<"KD22b: sum_triggeredAT_mirror:"<< sum_triggeredAT_mirror <<endl;
//++cout<<"KD22c: sum_triggeredST_mirror:"<< sum_triggeredST_mirror <<endl;

         }//end of mirror stations loop for st_type=3 or 4


         all_triggeredST = sum_triggeredST + sum_triggeredST_mirror;
         sum_ST = all_triggeredST;

//cout<<"KD  "<<"looks like I have finished surface and mirror calculations"<<endl;
//++cout<<"KD30: "<<"all_triggeredST= sum_triggeredST+sum_triggeredST_mirror= "<<all_triggeredST<<endl;
//++cout<<"------------------------------------------------------------------------------"<<endl;


         if (all_triggeredST < N_ST_required)
            continue;

         //cout<<"4.5th place"<<double(clock()-init_time)/CLOCKS_PER_SEC<<endl;
         count_events++;

         // cout<<"allTriST="<<all_triggeredST<<"  "<<sum_triggeredST<<"  "<<sum_triggeredST_mirror<<endl;
         int ff = 0;
         if (nuflavor == "nue")
            ff = 1;
         else if (nuflavor == "numu")
            ff = 2;
         else if (nuflavor == "nutau")
            ff = 3;
         else
            cout << "wrong nuflavor" << endl;

         int temp_totalAT = 0;
         int temp_totalAT_mirror = 0;
         int allAT = 0;


         if (sum_triggeredST > 0) {
            for (uint i = 0; i < N_TriggeredAT_perST_L2.size(); i++) { //loop over all the station triggered
               // cout<<N_TriggeredAT_perST_L2.at(i)<<endl;
               N_TriggeredAT_perST_L3.push_back(N_TriggeredAT_perST_L2.at(i));
               Trig_Type_L3.push_back(Trig_Type_L2.at(i));
               iRow_L3.push_back(iRow_L2.at(i));
               iCol_L3.push_back(iCol_L2.at(i));
               Efield_atposnu_L3.push_back(Efield_atposnu_L2.at(i));

               for (int j = 0; j < N_TriggeredAT_perST_L2.at(i); j++) //loop over all the triggered antennas in each station

               {
                  iAT_L3.push_back(iAT_L2.at(temp_totalAT + j));
                  path_inice_L3.push_back(path_inice_L2.at(temp_totalAT + j));
                  path_infirn_L3.push_back(path_infirn_L2.at(temp_totalAT + j));
                  totaltime_L3.push_back(totaltime_L2.at(temp_totalAT + j));
                  Efield_atposnu_L3.push_back(Efield_atposnu_L2.at(temp_totalAT + j));

                  volte_allAT_L3.push_back(volte_allAT_L2.at(temp_totalAT + j));
                  volth_allAT_L3.push_back(volth_allAT_L2.at(temp_totalAT + j));

                  AT_posi_x_L3.push_back(AT_posi_x_L2.at(temp_totalAT + j));

                  AT_posi_y_L3.push_back(AT_posi_y_L2.at(temp_totalAT + j));
                  AT_posi_z_L3.push_back(AT_posi_z_L2.at(temp_totalAT + j));

                  Dir_shower_x_L3.push_back(Dir_shower_x_L2.at(temp_totalAT + j));
                  Dir_shower_y_L3.push_back(Dir_shower_y_L2.at(temp_totalAT + j));
                  Dir_shower_z_L3.push_back(Dir_shower_z_L2.at(temp_totalAT + j));
                  theta_shower_atposnu_L3.push_back(theta_shower_atposnu_L2.at(temp_totalAT + j));

               }

               temp_totalAT += N_TriggeredAT_perST_L2.at(i);

            }


         }
         //cout<<"4.6th place"<<double(clock()-init_time)/CLOCKS_PER_SEC<<endl;
         if (sum_triggeredST_mirror > 0) {

            for (uint i = 0; i < N_TriggeredAT_perST_mirror_L2.size(); i++) { //loop over triggered stations

               N_TriggeredAT_perST_L3.push_back(N_TriggeredAT_perST_mirror_L2.at(i));
               Trig_Type_L3.push_back(Trig_Type_mirror_L2.at(i));
               iRow_L3.push_back(iRow_mirror_L2.at(i));
               iCol_L3.push_back(iCol_mirror_L2.at(i));



               for (int j = 0; j < N_TriggeredAT_perST_mirror_L2.at(i); j++) //loop over all the triggered antennas in each station

               {
                  // cout<<path_infirn_mirror_L2.at(temp_totalAT_mirror+j)<<endl;
                  iAT_L3.push_back(iAT_mirror_L2.at(temp_totalAT_mirror + j));
                  path_inice_L3.push_back(path_inice_mirror_L2.at(temp_totalAT_mirror + j));
                  path_infirn_L3.push_back(path_infirn_mirror_L2.at(temp_totalAT_mirror + j));
                  totaltime_L3.push_back(totaltime_mirror_L2.at(temp_totalAT_mirror + j));
                  Efield_atposnu_L3.push_back(Efield_atposnu_mirror_L2.at(temp_totalAT_mirror + j));


                  volte_allAT_L3.push_back(volte_allAT_mirror_L2.at(temp_totalAT_mirror + j));
                  volth_allAT_L3.push_back(volth_allAT_mirror_L2.at(temp_totalAT_mirror + j));

                  // cerr<<temp_totalAT_mirror+j<<endl;
                  //  cerr<<AT_posi_x_L2.at(temp_totalAT_mirror+j)<<endl;


                  AT_posi_x_L3.push_back(AT_posi_x_mirror_L2.at(temp_totalAT_mirror + j));
                  AT_posi_y_L3.push_back(AT_posi_y_mirror_L2.at(temp_totalAT_mirror + j));
                  AT_posi_z_L3.push_back(AT_posi_z_mirror_L2.at(temp_totalAT_mirror + j));

                  Dir_shower_x_L3.push_back(Dir_shower_x_mirror_L2.at(temp_totalAT_mirror + j));
                  Dir_shower_y_L3.push_back(Dir_shower_y_mirror_L2.at(temp_totalAT_mirror + j));
                  Dir_shower_z_L3.push_back(Dir_shower_z_mirror_L2.at(temp_totalAT_mirror + j));
                  theta_shower_atposnu_L3.push_back(theta_shower_atposnu_mirror_L2.at(temp_totalAT_mirror + j));

               }
               temp_totalAT_mirror += N_TriggeredAT_perST_mirror_L2.at(i);

            }

         }


//   //  cout<<"3335"<<endl;
//   // cout<<N_TriggeredAT_perST_L3.size()<<" stations triggered"<<endl;
         for (uint i = 0; i < N_TriggeredAT_perST_L3.size(); i++) //how many stations(including mirror stations) triggered in total=N_TriggeredAT_perST_L3.size();

         {
            // cout<<iRow_L3[i]<<"  "<<iCol_L3[i]<<endl;

            b3.NAT_atST[i] = N_TriggeredAT_perST_L3.at(i);
            b3.Trig_Type[i] = Trig_Type_L3.at(i);
            b3.ST_iRow[i] = iRow_L3.at(i);
            b3.ST_iCol[i] = iCol_L3.at(i);

            for (int j = 0; j < N_TriggeredAT_perST_L3.at(i); j++) {
               b3.iAT[allAT + j] = iAT_L3.at(allAT + j);
               b3.path_inice[allAT + j] = path_inice_L3.at(allAT + j);
               b3.path_infirn[allAT + j] = path_infirn_L3.at(allAT + j);
               b3.totaltime[allAT + j] = totaltime_L3.at(allAT + j);
               b3.Efield_atposnu[allAT + j] = Efield_atposnu_L3.at(allAT + j);
               b3.volte_allAT[allAT + j] = volte_allAT_L3.at(allAT + j);
               b3.volth_allAT[allAT + j] = volth_allAT_L3.at(allAT + j);
               b3.AT_Posi_x[allAT + j] = AT_posi_x_L3.at(allAT + j);
               b3.AT_Posi_y[allAT + j] = AT_posi_y_L3.at(allAT + j);
               b3.AT_Posi_z[allAT + j] = AT_posi_z_L3.at(allAT + j);
               b3.Dir_shower_x[allAT + j] = Dir_shower_x_L3.at(allAT + j);
               b3.Dir_shower_y[allAT + j] = Dir_shower_y_L3.at(allAT + j);
               b3.Dir_shower_z[allAT + j] = Dir_shower_z_L3.at(allAT + j);
               b3.theta_shower_atposnu[allAT + j] = theta_shower_atposnu_L3.at(allAT + j);

//      //cout<<b3.path_infirn[allAT+j]<<endl;
            }
            allAT += N_TriggeredAT_perST_L3.at(i);
         }
         b3.N_ATs = allAT;
//   // cout<<allAT<<" "<<iAT_L3.size()<<endl;//make sure they are equal

//      } //end of if(ST_TYPE=3 or 4)




//------KD added------7/8-------
      //for plotting
      theta_signal.push_back(theta_nposnu2AT * RAD2DEG); //theta of the signal when it first comes from the interaction point
      theta_signal_atAT.push_back(theta2 * RAD2DEG);
      viewangle_triggered.push_back(viewangle * RAD2DEG);
      attenfactor.push_back(exp(-d_posnu2AT / attenlength_up));
      dis.push_back(d_posnu2AT);
      vdis.push_back(ICETHICK - posnu[2]);
      cosz.push_back(cos(theta2));
      coszp.push_back(n_pol[2]);

      theta_signal_mirror.push_back((PI - theta_nposnu2MirrorAT)*RAD2DEG);
      theta_signal_atAT_mirror.push_back(theta2_mirror * RAD2DEG);
      viewangle_triggered_mirror.push_back(viewangle_mirror * RAD2DEG);
      attenfactor_mirror.push_back(exp(-d_posnu2MirrorAT / attenlength_down));
      dis_mirror.push_back(d_posnu2MirrorAT);
      vdis_mirror.push_back(ICETHICK + posnu[2]);
      cosz_mirror.push_back(cos(theta2_mirror));
      coszp_mirror.push_back(n_pol[2]);
//-----------------------------------------
//cout<<"dis.size()="<<dis.size()<<"  dis_mirror.size()="<<dis_mirror.size()<<endl;



      //cout<<"5th place"<<double(clock()-init_time)/CLOCKS_PER_SEC<<endl;

      weight = GetWeight(chord_inICE, chord_inCRUST, chord2_inICE, sigma, theta_nu);
//if (weight<0) {cout<<inu<<" "<<weight<<" ";}

//KD: 10/25 updating for tau regeneration
      if (TAUREGENERATION)      {
         if (nuflavor == "nutau" && log10(pnu) > 15. && log10(pnu) < 20. && theta_nu > (60 * DEG2RAD) && theta_nu < (PI / 2)) {
            weight = GetTauRegen(current, pnu, theta_nu, L_TauDecay, posnu[2]);
//      cout<<"ievt "<<inu<<" weight"<< weight<< endl;
         } else weight = GetWeight(chord_inICE, chord_inCRUST, chord2_inICE, sigma, theta_nu);

      }
//----------------------------------------------
//if (weight<0) //KD
//           {cout<<inu<<" "<<"tau weight "<<weight<<endl;}


      /*
      if (weight>0.01) //KD
      {
      outantposnu<< posnu[0] <<" "<< posnu[1] <<" "<<posnu[2]<<" "<<nnu[0] <<" "<< nnu[1] <<" "<<nnu[2]<< " " << weight <<" "<< ATCoordinate[0]<<" "<<ATCoordinate[1]<<" "<<ATCoordinate[2]<< " "<< abs_time << " "<< endl; //KD: it looks like it stops at NNU-1 and not for all NNU?
      }
      //outantposnu<< posnu2AT[0] <<" "<< posnu2AT[1] <<" "<<posnu2AT[2]<<" "<<nnu[0] <<" "<< nnu[1] <<" "<<nnu[2]<< " " << weight <<" "<< ATCoordinate[0]<<" "<<ATCoordinate[1]<<" "<<ATCoordinate[2]<< " "<< abs_time << " "<< endl; //KD: it looks like it stops at NNU-1 and not for all NNU?
      */

//if (weight>0.01) //KD
      integ_weight += weight;
//cout<<"integ_weight"<<integ_weight<<endl;//KD

      /*

      PUT BRACKETS HERE, LOOK AT ANDREA'S CODE, CHECK EVERYTHING

      */



//cout<<"viewangle_triggered.size()"<<viewangle_triggered.size()<<endl; //KD: getting zero here.

      b1.elpm     = elpm;
      b1.vmmhz1m_max = vmmhz1m_max;
      b1.changle  = changle * RAD2DEG;
      b1.eshower_em  = eshower_em;
      b1.eshower_had = eshower_had;

      b1.sum_triggeredST = sum_triggeredST;
      b1.sum_triggeredST_mirror = sum_triggeredST_mirror;

      b1.theta_nu = theta_nu * RAD2DEG;      //KD: hh1->Fill(theta_nu*RAD2DEG, weight);
      b1.depth = -(ICETHICK - posnu[2]);        //KD: depth->Fill(-(ICETHICK-posnu[2]),weight);
      b1.costhetanu = cos(theta_nu);         //KD: costhetanu->Fill(cos(theta_nu),weight);
      b1.phi_nu = phi_nu * RAD2DEG;

//      int ff = 0;
      if (nuflavor == "nue")
         ff = 1;
      else if (nuflavor == "numu")
         ff = 2;
      else if (nuflavor == "nutau")
         ff = 3;

      int mycurrent = 0;
      if (current == "cc")
         mycurrent = 1;
      else if (current == "nc")
         mycurrent = 2;


      b1.ievt = inu;
      b1.Energy = EXPONENT;
      b1.flavor = ff;
      b1.mycurrent = mycurrent;
      b1.y = elast_y;
      b1.weight = weight;
      b1.SumHits = sum_ST; //number of stations with greater than 0.5NV signal //KD: not valid for ST_TYPE=4 or 3?

//-------------------------------
      b7.N_Ant_perST = N_Ant_perST;
//-------------------------------

      //b1.N_allchannels=n_chan_perST*sum_ST; //KD: valid for ST_TYPE 0,1,2
      b1.N_Ant_perST = N_Ant_perST;
      b1.NTriggeredST = all_triggeredST;
      b1.Dir_nu[0] = nnu[0];
      b1.Dir_nu[1] = nnu[1];
      b1.Dir_nu[2] = nnu[2];
      b1.Posi_Int[0] = posnu[0];
      b1.Posi_Int[1] = posnu[1];
      b1.Posi_Int[2] = posnu[2];

      //KD: added 01/13/11, rudimentary for 1 station, needs improvement
      if (PLANEWAVE) {
         //for(uint i=0; i<sum_triggeredST;i++){
         if (sum_triggeredST == 1) {
            b1.my_ST_Posi[0][0] = ATCoordinate[0];
            b1.my_ST_Posi[0][1] = ATCoordinate[1];
            b1.my_ST_Posi[0][2] = ATCoordinate[2];
            b1.my_d_posnu2ST = d_posnu2ST;
            b1.my_ST_abs_time = ST_abs_time;
            b1.my_nsignal_atST[0][0] = nsignal_atST[0];
            b1.my_nsignal_atST[0][1] = nsignal_atST[1];
            b1.my_nsignal_atST[0][2] = nsignal_atST[2];

         }

         if (sum_triggeredST_mirror == 1) {
            b1.my_ST_Posi[1][0] = MirrorATCoordinate[0];
            b1.my_ST_Posi[1][1] = MirrorATCoordinate[1];
            b1.my_ST_Posi[1][2] = MirrorATCoordinate[2];
            b1.my_d_posnu2MirrorST = d_posnu2MirrorST;
            b1.my_ST_abs_time_mirror = ST_abs_time_mirror;
            b1.my_nsignal_atST[1][0] = nsignal_atMirrorST[0];
            b1.my_nsignal_atST[1][1] = nsignal_atMirrorST[1];
            b1.my_nsignal_atST[1][2] = nsignal_atMirrorST[2];
         }

      }


      //KD: me filling new ones
      b1.viewangle_triggered = viewangle * RAD2DEG;      //vangle->Fill(viewangle_triggered.at(i),weight);
      b1.viewangle_triggered_mirror = viewangle_mirror * RAD2DEG; //vangle_mirror->Fill(viewangle_triggered_mirror.at(i));
      b1.attenfactor = exp(-d_posnu2AT / attenlength_up);      //dattenfactor->Fill(attenfactor.at(i),weight);
      b1.attenfactor_mirror = exp(-d_posnu2MirrorAT / attenlength_down); //rattenfactor->Fill(attenfactor_mirror.at(i),weight);
      b1.dis = d_posnu2AT;                //dir_dis->Fill(dis.at(i),weight);
      b1.dis_mirror = d_posnu2MirrorAT;               //ref_dis->Fill(dis_mirror.at(i),weight);
      b1.vdis = ICETHICK - posnu[2];            //dir_vdis->Fill(vdis.at(i),weight);
      b1.vdis_mirror = ICETHICK + posnu[2];        //ref_vdis->Fill(vdis_mirror.at(i),weight);
      b1.cosz = cos(theta2);                 //dir_cosz->Fill(cosz.at(i),weight);
      b1.cosz_mirror = cos(theta2_mirror);            //ref_cosz->Fill(cosz_mirror.at(i),weight);
      b1.coszp = n_pol[2];    //KD 11/15: not correct? calls the last mirror value                 //dir_coszp->Fill(coszp.at(i),weight);
      b1.coszp_mirror = n_pol[2];//KD 11/15: not correct? calls the last mirror value              //ref_coszp->Fill(coszp_mirror.at(i),weight);

      //KD: added on 01/11/11
      b1.theta_signal = theta_nposnu2AT * RAD2DEG;
      //b1.theta_signal_atAT = theta2*RAD2DEG;
      b1.theta_signal_mirror = (PI - theta_nposnu2MirrorAT) * RAD2DEG;
      // b1.theta_signal_atAT_mirror = theta2_mirror*RAD2DEG;

      //b1.phi_nposnu2AT = phi_nposnu2AT*RAD2DEG;
      //b1.phi_nposnu2MirrorAT = phi_nposnu2MirrorAT*RAD2DEG;


// KD 11/15/2010     !!!!!!!! wrong to put it here, as it calls the mirror values last.
//     b1.Polarization[0] = n_pol[0];
//     b1.Polarization[1] = n_pol[1];
//     b1.Polarization[2] = n_pol[2]; //KD: redundant with coszp and coszp_mirror above

//      SphAngles(n_pol, sphangles);
//      b1.polariz_phi = sphangles[0] * RAD2DEG;
// b1.polariz_theta = sphangles[1] * RAD2DEG;


      GetPolarization(nnu, nposnu2AT, n_pol, n_Bfield);
      b1.Polariz_direct[0] = n_pol[0];
      b1.Polariz_direct[1] = n_pol[1];
      b1.Polariz_direct[2] = n_pol[2];
      b1.Bfield_direct[0] = n_Bfield[0];
      b1.Bfield_direct[1] = n_Bfield[1];
      b1.Bfield_direct[2] = n_Bfield[2];


      SphAngles(n_pol, sphangles);
      b1.Polariz_direct_phi = sphangles[0] * RAD2DEG;
      b1.Polariz_direct_theta = sphangles[1] * RAD2DEG;


      GetPolarization(nnu, nposnu2MirrorAT, n_pol, n_Bfield);
      b1.Polariz_reflect[0] = n_pol[0];
      b1.Polariz_reflect[1] = n_pol[1];
      b1.Polariz_reflect[2] = n_pol[2];
      b1.Bfield_reflect[0] = n_Bfield[0];
      b1.Bfield_reflect[1] = n_Bfield[1];
      b1.Bfield_reflect[2] = n_Bfield[2];


      SphAngles(n_pol, sphangles);
      b1.Polariz_reflect_phi = sphangles[0] * RAD2DEG;
      b1.Polariz_reflect_theta = sphangles[1] * RAD2DEG;


      tree1->Fill();
//-------------------------------
      tree7->Fill();
//-------------------------------

//------------------added KD for multi-st-----------------
      b3.ievt = inu;
      b3.Energy = EXPONENT;
      b3.flavor = ff;
      b3.y = elast_y;
      b3.weight = weight;
      b3.NTriggeredST = all_triggeredST;
      b3.Theta_nu = theta_nu;
      b3.Dir_nu[0] = nnu[0];
      b3.Dir_nu[1] = nnu[1];
      b3.Dir_nu[2] = nnu[2];
      b3.Posi_Int[0] = posnu[0];
      b3.Posi_Int[1] = posnu[1];
      b3.Posi_Int[2] = posnu[2];

      b3.sum_triggeredST = sum_triggeredST;
      b3.sum_triggeredST_mirror = sum_triggeredST_mirror;

      tree3->Fill();
//      if(weight>0.001)


//------------------added KD for multi-st-----------------

// separating neutrino flavor for Veff studies
      if (nuflavor == "nue") {
         events_nue++; //KD
         nue_w += weight; //KD
      } else if (nuflavor == "numu") {
         events_numu++; //KD
         numu_w += weight; //KD
      } else {
         events_nutau++; //KD
         nutau_w += weight; //KD
      }

//-----------------------------------------------------------------------

      // different definitions of  direct events or reflected events,or combo events for fenfang's plots
      if (sum_triggeredST == 0) {
         //fan_ref->Fill(EXPONENT,weight);
         reflect_events++; //KD
         fan_ref_w += weight; //KD
         //thr->Fill(theta_nu,weight); //KD
      } else if (sum_triggeredST_mirror == 0) {
         //fan_dir->Fill(EXPONENT,weight);
         direct_events++; //KD
         fan_dir_w += weight; //KD
         //thd->Fill(theta_nu*RAD2DEG,weight); //KD
         //thd->Fill(theta_nu*RAD2DEG); //KD

      } else {
//KD: this seems to get populated if both direct and reflected hits the station
         //fan_combo->Fill(EXPONENT,weight);
         combo_events++; //KD
         fan_combo_w += weight; //KD
         // thc->Fill(theta_nu,weight); //KD
      }

      //  tree0->Fill();

   }//end of neutrino loop //KD: YEP THIS IS THE END OF THE NNU LOOP
   fprintf(stderr,"\n");
//===============================================================

   // cout<<dir->GetBin(17)<<endl;

   // cout<<dir->GetBinContent(0)<<"  "<<dir->GetBinContent(1)<<" "<<dir->GetBinContent(7)<<endl;

   // cout<<sum_bin_content<<"  "<<integ_weight<<endl;
   // output<<"count_events"<<count_events<<endl;
   // Veff=volume*DensityICE/DensityWATER*4*PI/NNU*integ_weight/2*1.E-9;//km3sr
   Veff = volume * DensityICE / DensityWATER * 4 * PI / NNU * integ_weight * 1.E-9;
   //km3sr, start to use it since 6:00pm July 13th, 2006

//KD reactivated line below on Nov 2011, looks like livetime is already in seconds?
//KD I'm not sure equation is used in the right way here
   E2dNdE = 1.e-4 * 2.3 * int_len_kgm2 / DensityWATER / (Veff * 1.e9) / LIVETIME * 1.e10 * YR2SEC; //the unit is cm^-2 sr^-1 s^-1

   Veff_nue = volume * DensityICE / DensityWATER * 4 * PI / nue_counts * nue_w * 1.E-9;
   Veff_numu = volume * DensityICE / DensityWATER * 4 * PI / numu_counts * numu_w * 1.E-9;
   Veff_nutau = volume * DensityICE / DensityWATER * 4 * PI / nutau_counts * nutau_w * 1.E-9;
   //cout<<endl;
   //cout<<"==================================="<<endl;

   cout << "NV= " << NV << " REFLECT_RATE=" << REFLECT_RATE << " ATTEN_UP=" << ATTEN_UP << " ICETHICK=" << ICETHICK << " FIRNfactor=" << FIRNfactor << endl;

   output << "NV= " << NV << " REFLECT_RATE=" << REFLECT_RATE << " ATTEN_UP=" << ATTEN_UP << " ICETHICK=" << ICETHICK << " FIRNfactor=" << FIRNfactor << endl;


   if (PLANEWAVE) {
      //kd added for detailed output
      double ATCoordinates8[N_Ant_perST][3];
      double MirrorATCoordinates8[N_Ant_perST][3];

      for (int i = 0; i < N_Ant_perST; i++) {
         double phi = (2. / N_Ant_perST) * PI * i; //the phi angle of each LPA's center
         ATCoordinates8[i][0] = ATCoordinate[0] + ST4_R * cos(phi);
         ATCoordinates8[i][1] = ATCoordinate[1] + ST4_R * sin(phi);
         ATCoordinates8[i][2] = ATCoordinate[2];
//cout<<"KD5i: "<<"ATCoordinates: ("<<ATCoordinates8[i][0]<<","<<ATCoordinates8[i][1]<<","<<ATCoordinates8[i][2]<<")"<<endl;

         //phi=(2./N_Ant_perST)*PI*i;//the phi angle of each LPA's center
         MirrorATCoordinates8[i][0] = MirrorATCoordinate[0] + ST4_R * cos(phi);
         MirrorATCoordinates8[i][1] = MirrorATCoordinate[1] + ST4_R * sin(phi);
         MirrorATCoordinates8[i][2] = MirrorATCoordinate[2];
//cout<<"KD5ii: "<<"MirrorATCoordinates: ("<<MirrorATCoordinates8[i][0]<<","<<MirrorATCoordinates8[i][1]<<","<<MirrorATCoordinates8[i][2]<<")"<<endl;
      }
   }


   cout <<
        "reflect_events,weighted= " << reflect_events << "," << fan_ref_w <<
        " direct_events,weighted= " << direct_events << "," << fan_dir_w <<
        " combo_events,weighted= " << combo_events << "," << fan_combo_w <<
        endl;

   cout <<
        "nue/counts,weighted, Veff_nue= " << events_nue << "/" << nue_counts << "," << nue_w << "," << Veff_nue <<
        " numu/counts,weighted, Veff_numu= " << events_numu << "/" << numu_counts << "," << numu_w << "," << Veff_numu <<
        " nutau/counts,weighted, Veff_nutau= " << events_nutau << "/" << nutau_counts << "," << nutau_w << "," << Veff_nutau <<
        endl;



   output <<
        "reflect_events,weighted= " << reflect_events << "," << fan_ref_w <<
        " direct_events,weighted= " << direct_events << "," << fan_dir_w <<
        " combo_events,weighted= " << combo_events << "," << fan_combo_w <<
        endl;

   output <<
        "nue/counts,weighted, Veff_nue= " << events_nue << "/" << nue_counts << "," << nue_w << "," << Veff_nue <<
        " numu/counts,weighted, Veff_numu= " << events_numu << "/" << numu_counts << "," << numu_w << "," << Veff_numu <<
        " nutau/counts,weighted, Veff_nutau= " << events_nutau << "/" << nutau_counts << "," << nutau_w << "," << Veff_nutau <<
        endl;


   cout << "Veff=" << Veff << " km3sr; integ_weight=" << integ_weight << " NNU=" << NNU << " volume=" << volume << " seed=" << seed << "  AT/ST= " << N_Ant_Trigger << "/" << N_Ant_perST << "  NSIGMA=" << NSIGMA << " gainv=" << gainv << " " << "array=" << NROWS << "x" << NCOLS << " E2dNdE=" << E2dNdE << " ST4_R=" << ST4_R << endl;

//   output << "Veff= " << Veff << " km3sr  " << "NNU= " << NNU << " volume=" << volume << " seed= " << seed << " AT/ST= " << N_Ant_Trigger << "/" << N_Ant_perST << " integ_weight= " << integ_weight << " " << "  NSIGMA=" << NSIGMA << " gainv=" << gainv << " " << " array=" << NROWS << "x" << NCOLS << " R=" << REFLECT_RATE << " ATTEN_UP=" << ATTEN_UP;

   output << "Veff=" << Veff << " km3sr; integ_weight=" 
          << integ_weight << " NNU=" << NNU << " volume=" << volume 
          << " seed=" << seed << "  AT/ST= " << N_Ant_Trigger << "/" 
          << N_Ant_perST << "  NSIGMA=" << NSIGMA << " gainv=" 
          << gainv << " " << "array=" << NROWS << "x" << NCOLS 
          << " E2dNdE=" << E2dNdE << " ST4_R=" << ST4_R 
          << " R=" << REFLECT_RATE 
          <<" ATTEN_UP=" << ATTEN_UP<< endl;

//" ICETHICK="<<ICETHICK;
// KD no 'endl' for concise output
//output<<"Veff= "<<Veff<<" km3sr  "<<endl;

   if (!SPECTRUM) {
      cout << pnu << "  " << Veff << " km3sr" << endl;
      output << " pnu= " << pnu << "  " << Veff << " km3sr" << endl;
//    outveff<<pnu<<"  "<<Veff<<" km3sr"<<endl; //KD: I deactivated outveff file
   }

   CloseTFile(hfile);


   if (SPECTRUM)

   {
      double sum_events = 0;
      double thisenergy = 0;
      double thislen_int_kgm2 = 0;
      if (WIDESPECTRUM) {
         for (int i = 0; i < 12; i++) {
            EdNdEdAdt[i] = EdNdEdAdt[i] * maxflux; //maxflux gotten the main inu loop, it's the same for all
         } //for
         for (int i = 0; i < 12; i++) {
            thisenergy = pow(10., 16. + ((double)i) * 0.5);
            Sigma(thisenergy, sigma, thislen_int_kgm2);
            sum_events += EdNdEdAdt[i] / (thislen_int_kgm2 / DensityWATER);
         } //for

      } else {
         for (int i = 0; i < 7; i++) {
            EdNdEdAdt[i] = EdNdEdAdt[i] * maxflux; //maxflux gotten the main inu loop, it's the same for all
         } //for
         for (int i = 0; i < 7; i++) {
            thisenergy = pow(10., 17. + ((double)i) * 0.5);
            Sigma(thisenergy, sigma, thislen_int_kgm2);
            sum_events += EdNdEdAdt[i] / (thislen_int_kgm2 / DensityWATER);
         } //for
      }
      // sum_events*=volume*LIVETIME*RHOICE/RHOH20*nevents/NNU*sr*pow(10.,4);
      sum_events *= Veff * LIVETIME * pow(10, 13.) * log(10.) * 0.5;
      //  cout<<volume*LIVETIME*DensityICE/DensityWATER*integ_weight/NNU*4.*PI*pow(10,4.)<<"    =?";
      //KD 11/9/11: 0.5 is the bin size of integration d(log E)???

      cout << sum_events;

      if (FANFLUX)
         sum_events *= 0.75; //using E-2 flux that I was using always but times 0.75
      // cout<<Veff*LIVETIME*pow(10,13.)<<endl;
      cout << "   sum of events is " << sum_events << " pnu,maxflux:" << pnu << "," << maxflux << "\n";

      output << " sum of events= " << sum_events << endl; //KD for easy extraction for stats test
      output << "   sum of events is " << sum_events << " pnu,maxflux:" << pnu << "," << maxflux << "\n";

      //   output<<"ATGap"<<"  "<<"ST_TYPE"<<"  "<<"FIRN"<<"  "<<"NFIRN"<<"  "<<"SCATTER"<<"  "<<"SCATTER_WIDTH"<<"  "<<"CONST_ATTENLENGTH"<<"  "<<"NSIGMA"<<"  "<<"sum_events"<<endl;

//FW      output<<ATGap<<"  "<<ST_TYPE<<"  "<<FIRN<<"  "<<NFIRN<<"  "<<SCATTER<<"  "<<SCATTER_WIDTH<<"  "<<CONST_ATTENLENGTH<<"  "<<NSIGMA<<"  "<<sum_events<<endl;

      // output<<JAIME_FACTOR<<"   "<<sum_events<<endl;
      //  outveff<<ATTEN_UP<<"  "<<REFLECT_RATE<<" "<<sum_events<<endl;
//     outveff<<ATGap<<"  "<<sum_events<<endl;

   }//end of if(SPECTRUM)


//KD: End of output files






   //KD: END OF OUTPUT FILES



//--DISPLAYING END TIME OF PROGRAM------------------------+
   time_t raw_final_time = time(NULL);       //|
   struct tm* final_time = localtime(&raw_final_time);   //|
//clock_t final_time=clock();          //|
   cout << asctime(final_time) << endl;   //|
//--------------------------------------------------------+

   return 0;

}

