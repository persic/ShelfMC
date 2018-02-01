#ifndef __antennamodel_ara_h__   // first line in file
#define __antennamodel_ara_h__   // second line in file

#include <iostream>
#include <string>
#include <TChain.h>
#include <TArrayF.h>
#include <TVector3.h>

using namespace std;

class ARA_Ant;
class TChain;

class ARA_Ant {

    private:
    
        TChain* nt;
        TChain* ntt;
    
    public:
    
        ARA_Ant();
        char* filename;
//         int nfreq; // Number of frequencies that tree was made with
        ARA_Ant(char* filename);
        
        // AntTree
        Int_t N[1];
        Float_t thetas[1];
        Float_t phis[1];
        Float_t frequencies[60];
        Double_t gains[60];
        Double_t phases[60];

        
        // ZTree
        Int_t NN[0];
        Double_t Re_Z[60];
        
        double theta_d;
        double phi_d;
        
        int theta_i;
        int phi_i;
        int n;
        
        void LoadData(double theta_d, double phi_d);
        void CheckBoundaries(double &theta_d, double &phi_d);
        TVector3* GetVector(double zenith,double azimuth); 
        
        struct AntennaAngles {
            double theta;
            double phi;
            };
            
//         AntennaAngles GetThetaAndPhiFlexible(double zenith, 
//                                              double azimuth, 
//                                              double azi_orientation, 
//                                              double zen_orientation);
        
        double InterpolateToSingleFrequency(double freq, int nfreq, float* frequencies, double* gains);
//         
        void LoadGain(double* n_dipole, 
                        double* n_arrivaldir); 
//                             
        double GetEffectiveHeight( 
                                double gain,
                                double freq,
                                double ant_Z, 
                                double c,
                                double Z_0);
        
        virtual ~ARA_Ant(); //destructor
};

#endif  // last line in file
