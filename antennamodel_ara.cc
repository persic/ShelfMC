#include <string>
#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TArray.h>
#include <TVector3.h>
#include <TMath.h>

#include "antennamodel_ara.h"

ARA_Ant::ARA_Ant()
    {
    //Default constructor
    }
    
    
ARA_Ant::ARA_Ant(char* filename)

{
    nt = new TChain("AntTree");
    nt->Add(filename);
    
    nt->SetBranchAddress("N", &N);
    nt->SetBranchAddress("thetas", &thetas);
    nt->SetBranchAddress("phis", &phis);
    nt->SetBranchAddress("frequencies", &frequencies);
    nt->SetBranchAddress("gains", &gains);
    nt->SetBranchAddress("phases", &phases);

    nt->BuildIndex("thetas", "phis"); 

    ntt = new TChain("ZTree");
    ntt->Add(filename);
    
    ntt->SetBranchAddress("N", &NN);        
    ntt->SetBranchAddress("Re_Z", &Re_Z);
    }

void ARA_Ant::LoadData(double theta_d, double phi_d){

    int theta = static_cast<int>(theta_d);
    int phi = static_cast<int>(phi_d);
    
    int mod_theta = theta % 5;
    int mod_phi = phi % 5;
    
    // ARA antenna model only provides angles in steps of 5 degrees.
    if (mod_theta <= 2){
        theta -= mod_theta;
        }
    else {
        theta += (5 - mod_theta);
    }
    
    if (mod_phi <= 2){
        phi -= mod_phi;
        }
    else {
        phi += (5 - mod_phi);
    }  
    
    if (phi == 360){
        phi = 0;
        }
    
    // Load impedance tree
    int m = ntt->GetEntry(0);
    if (m == -1){
        cout << "Impedance tree not filled" << endl;
        }
        
    // Load Gain tree
    int n = nt->GetEntryNumberWithIndex(theta,phi);
    if (n > -1){
        nt->GetEntry(n);
        for (int i = 0; i < N[0]; i++){
            // Convert frequencies to MHz as default
            frequencies[i] = frequencies[i]* 1000;
            }
        }    
    else {
        cout << "Direction theta "<< theta << " and phi " << phi << " not found in Tree" << endl;
        exit(1);
        }     

    }
    
void ARA_Ant::CheckBoundaries(double &theta_d, double &phi_d){
    // Wrapping of angles according to angel convention
    if (theta_d > 180){
        theta_d -= 180;
        }
    if (theta_d < 0){
        theta_d += 180;
        }
    if (phi_d < 0){
        phi_d += 360;
        }
    if (phi_d >= 360){
        phi_d -= 360;
        }
    }
    
TVector3*  ARA_Ant::GetVector(double zenith,double azimuth){
    //Initialize vector for rotation.
    TVector3* v = new TVector3(1,0,0);
    v->SetPhi(azimuth/180. * TMath::Pi() );
    v->SetTheta(zenith/180. * TMath::Pi());

    return v;  
    }
 
 
// ARA_Ant::AntennaAngles ARA_Ant::GetThetaAndPhiFlexible(double zenith, 
//                                                  double azimuth, 
//                                                  double azi_orientation=90., 
//                                                  double zen_orientation=90.){
        /// Get antenna angles for various rotations, azimuth along tines"""
//         TVector3* v = GetVector(zenith, azimuth);
//         v->RotateZ(-1*azi_orientation * TMath::Pi() / 180.);
//         v->RotateY((90-zen_orientation) * TMath::Pi() / 180.);
//         AntennaAngles ant_angles = Cartesian2WPLD(v);
//         return ant_angles;
//         }

double ARA_Ant::InterpolateToSingleFrequency(double freq, int nfreq, float* frequencies, double* gains){
        // provide interpolated value of gains for single frequency freq in nfreq frequencies
        int i = 0;
        for (int j = 0; j < nfreq; j++ ){
            if (frequencies[j] > freq) {
                i = j;
                break;
                }
            }
        double gain = 0;
        if (frequencies[i] == freq){
            gain = gains[i];
            }
        else {
            // simple linear interpolation
            gain = (gains[i] - gains[i-1])/(frequencies[i] - frequencies[i-1]) * (freq - frequencies[i-1]) + gains[i-1];
            }
       
        return gain;
        }


double ARA_Ant::GetEffectiveHeight(double gain,
                                double freq,
                                double ant_Z, 
                                double c=2.99792e8,
                                double Z_0=119.9169*TMath::Pi()){

            // Antenna gain converted to Antenna effective height, see also thesis Kamlesh """
            double A_e = TMath::Sqrt(gain/TMath::Pi()*TMath::Power(c/(freq*1e6),2)*(ant_Z/Z_0));
            
        return A_e;
    } 
        
void ARA_Ant::LoadGain( double* n_dipole, 
                        double* n_arrivaldir){
                                  
      TVector3 dipole(n_dipole[0],n_dipole[1],n_dipole[2]); 
      TVector3 arrivaldir(n_arrivaldir[0],n_arrivaldir[1],n_arrivaldir[2]); 
                  
      double zen_orientation = dipole.Theta()/ TMath::Pi() * 180.;
      double azi_orientation = dipole.Phi()/ TMath::Pi() * 180.;
      
      double zenith = arrivaldir.Theta()/ TMath::Pi() * 180.;
      double azimuth = arrivaldir.Phi()/ TMath::Pi() * 180.;
            
      TVector3* v = GetVector(zenith, azimuth);
        
      v->RotateZ(-1*azi_orientation * TMath::Pi() / 180.);
      v->RotateY((-1*zen_orientation) * TMath::Pi() / 180.);
           
      AntennaAngles ant_angles;
      ant_angles.theta = v->Theta()/ TMath::Pi() * 180.;
      ant_angles.phi = v->Phi()/ TMath::Pi() * 180.;
      
//       cout << "1) Angles to be accessed, theta " << ant_angles.theta << " phi " << ant_angles.phi << endl;
      
      CheckBoundaries(ant_angles.theta, ant_angles.phi);
      
//       cout << "2) Angles to be accessed, theta " << ant_angles.theta << " phi " << ant_angles.phi << endl;
      LoadData(ant_angles.theta,ant_angles.phi);
                                                             
     }


 
ARA_Ant::~ARA_Ant() {
    delete nt;
    delete ntt;
}   
