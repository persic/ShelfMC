#ifndef __antennamodel_lpdaZ_h__   // first line in file
#define __antennamodel_lpdaZ_h__   // second line in file

#include <iostream>
#include <string>
#include <TChain.h>
#include <TArrayF.h>
#include <TVector3.h>

using namespace std;

class LPDAZ;
class TChain;

class LPDAZ {

    private:

        TChain* nt;
        TChain* ntt;

    public:

        LPDAZ();
        char* filename;
//         int nfreq; // Number of frequencies that tree was made with
        LPDAZ(char* filename);

        // AntTree
        Int_t N[1];
        Float_t thetas[1];
        Float_t phis[1];
        Float_t frequencies[219];
        Double_t gains[219];
        Double_t Re_phi[219];
        Double_t Im_phi[219];
        Double_t Re_theta[219];
        Double_t Im_theta[219];

        // ZTree
        Int_t NN[0];
        Double_t Re_Z[219];
        Double_t Im_Z[219];

        double theta_d;
        double phi_d;

        int theta_i;
        int phi_i;
        int n;

        void LoadData(double theta_d, double phi_d);
        void CheckBoundaries(double &theta_d, double &phi_d);
        TVector3* GetVector(double zenith,double azimuth);
        TVector3* GetPhiHat(TVector3* v);
        TVector3* GetThetaHat(TVector3* v);

        struct AntennaAngles {
            double theta;
            double phi;
            };

        AntennaAngles Cartesian2WPLD(TVector3* v);
        AntennaAngles WPLD2Cartesian(TVector3* v);
        AntennaAngles GetThetaAndPhiFlexible(double zenith,
                                             double azimuth,
                                             double azi_orientation,
                                             double zen_orientation);

        double InterpolateToSingleFrequency(double freq, int nfreq, float* frequencies, double* gains);

        void LoadGain(double* n_boresight,
                        double* n_eplane,
                        double* n_arrivaldir);

        double GetEffectiveHeight(
                                double gain,
                                double freq,
                                double ant_Z,
                                double c,
                                double Z_0);

        double GetEffectiveLength(double Z,
                                        double wavelength,
                                        double I_phi,
                                        double I_theta,
                                        double* n_boresight,
                                        double* n_eplane,
                                        double* n_arrivaldir,
                                        double* n_pol,
                                        double Z_0);

        TVector3* AlignVector(double* n_boresight,
                                  double* n_eplane,
                                  double* n_arrivaldir);

        virtual ~LPDAZ(); //destructor
};

#endif  // last line in file