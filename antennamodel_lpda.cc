#include <string>
#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TArray.h>
#include <TVector3.h>
#include <TMath.h>

#include "antennamodel_lpda.h"

LPDA::LPDA()
    {
    //Default constructor
    }
    
    
LPDA::LPDA(char* filename)

{
    nt = new TChain("AntTree");
    nt->Add(filename);
    
    nt->SetBranchAddress("N", &N);
    nt->SetBranchAddress("thetas", &thetas);
    nt->SetBranchAddress("phis", &phis);
    nt->SetBranchAddress("frequencies", &frequencies);
    nt->SetBranchAddress("gains", &gains);
    //    nt->SetBranchAddress("Re_phi", &Re_phi);
    //    nt->SetBranchAddress("Im_phi", &Im_phi);
    //    nt->SetBranchAddress("Re_theta", &Re_theta);
    //    nt->SetBranchAddress("Im_theta", &Im_theta);

    nt->BuildIndex("thetas", "phis"); 

    ntt = new TChain("ZTree");
    ntt->Add(filename);
    
    ntt->SetBranchAddress("N", &NN);        
    ntt->SetBranchAddress("Re_Z", &Re_Z);
    ntt->SetBranchAddress("Im_Z", &Im_Z);
    }

void LPDA::LoadData(double theta_d, double phi_d){

    int theta = static_cast<int>(theta_d);
    int phi = static_cast<int>(phi_d);
    
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
        cout << "Direction theta "<< theta << " and phi " << phi << " not found in WiplD" << endl;
        exit(1);
        }     

    }
    
void LPDA::CheckBoundaries(double &theta_d, double &phi_d){
    // Wrapping of angles according to angel convention of WIPL-D.
    if (theta_d > 90){
        theta_d -= 180;
        }
    if (theta_d < -90){
        theta_d += 180;
        }
    if (phi_d < 0){
        phi_d += 360;
        }
    if (phi_d > 360){
        phi_d -= 360;
        }
    }
    
TVector3*  LPDA::GetVector(double zenith,double azimuth){
    //Initialize vector for rotation.
    TVector3* v = new TVector3(1,0,0);
    v->SetPhi(azimuth/180. * TMath::Pi() );
    v->SetTheta(zenith/180. * TMath::Pi());

    return v;  
    }

LPDA::AntennaAngles LPDA::Cartesian2WPLD(TVector3* v){
        /// Conversion from Cartesian to WiPLD
        double theta = 90. - v->Theta()/TMath::Pi() * 180.;
        double phi = v->Phi()/TMath::Pi() * 180.;
        
        CheckBoundaries(theta, phi);
        AntennaAngles ant_angles;
        ant_angles.theta = theta;
        ant_angles.phi = phi;
        return ant_angles;
        }   

LPDA::AntennaAngles LPDA::WPLD2Cartesian(TVector3* v){
        /// Conversion from WiPLD to Cartesian (same conversion as other)
        AntennaAngles ant_angles = Cartesian2WPLD(v);
        return ant_angles;
        } 
 
LPDA::AntennaAngles LPDA::GetThetaAndPhiFlexible(double zenith, 
                                                 double azimuth, 
                                                 double azi_orientation=90., 
                                                 double zen_orientation=90.){
        /// Get antenna angles for various rotations, azimuth along tines"""
        TVector3* v = GetVector(zenith, azimuth);
        v->RotateZ(-1*azi_orientation * TMath::Pi() / 180.);
        v->RotateY((90-zen_orientation) * TMath::Pi() / 180.);
        AntennaAngles ant_angles = Cartesian2WPLD(v);
        return ant_angles;
        }

double LPDA::InterpolateToSingleFrequency(double freq, int nfreq, float* frequencies, double* gains){
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


double LPDA::GetEffectiveHeight(double gain,
                                double freq,
                                double ant_Z, 
                                double c=2.99792e8,
                                double Z_0=119.9169*TMath::Pi()){

            // Antenna gain converted to Antenna effective height, see also thesis Kamlesh """
            double A_e = TMath::Sqrt(gain/TMath::Pi()*TMath::Power(c/(freq*1e6),2)*(ant_Z/Z_0));
            
        return A_e;
    } 
        
void LPDA::LoadGain(  double* n_boresight, 
                      double* n_eplane, 
                      double* n_arrivaldir){
                                  
      TVector3 boresight(n_boresight[0],n_boresight[1],n_boresight[2]); 
      TVector3 eplane(n_eplane[0],n_eplane[1],n_eplane[2]);  
      TVector3 arrivaldir(n_arrivaldir[0],n_arrivaldir[1],n_arrivaldir[2]);                            
      TVector3 hplane = eplane.Cross(boresight);
                                  
      double zen_orientation = boresight.Theta()/ TMath::Pi() * 180.;
      double azi_orientation = boresight.Phi()/ TMath::Pi() * 180.;
      
      double zenith = arrivaldir.Theta()/ TMath::Pi() * 180.;
      double azimuth = arrivaldir.Phi()/ TMath::Pi() * 180.;
      
      eplane.RotateZ(-1*azi_orientation * TMath::Pi() / 180.);
      eplane.RotateY((90-zen_orientation) * TMath::Pi() / 180.);
      double additional_rotation = eplane.Theta();
            
      TVector3* v = GetVector(zenith, azimuth);
        
      v->RotateZ(-1*azi_orientation * TMath::Pi() / 180.);
      v->RotateY((90-zen_orientation) * TMath::Pi() / 180.);
      v->RotateX(additional_rotation-TMath::Pi()/2.);
           
      AntennaAngles ant_angles = Cartesian2WPLD(v);
      
      CheckBoundaries(ant_angles.theta, ant_angles.phi);

      LoadData(ant_angles.theta,ant_angles.phi);
                                                             
     }


 
LPDA::~LPDA() {
    delete nt;
    delete ntt;
}   
