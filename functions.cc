
void GetNuFlavor(string& nuflavor)
{
   double temp = Rand3.Rndm();

   if (temp <= (1. / 3.)) {
      nuflavor = "nue";
   } else if (temp <= (2. / 3.)) {
      nuflavor = "numu";
   } else {
      nuflavor = "nutau";

   }

}

void GetCurrent(string& current)
{
   //CC or NC
   double temp = Rand3.Rndm();
   if (temp <= 0.7064)
      current = "cc";
   else
      current = "nc";


}


/*int RossIceshelf(double *position)
  {
  int lon,lat;
  GetLonLat(postion,lon,lat);//West Longitude is negative
  if(((lon<=180&&lon>=150&&lat>=78&&lat<=86)!
  (lat>=78&&lat<=84&&lon<=160&&lon>=150))||
  (lon>=-180&&lon<=-150&&lat>=78&&lat<=84)||
  (lon>=-180&&lon<=-170&&lat>=84&&lat<=86)||
  (lon>=-150&&lon<=-140&&lat>=80&&lat<=82)
  )
  return 1;
  else return 0;

  }*/


//set the junction of E180 longitude and S80 latitude as the origin
// E180 northward as y, S80 westward as x, down to the ice as z

void  GetBothLocation(int iRow, int iCol, double* ATCoordinate, double* MirrorATCoordinate)
{
   if (HEXAGONAL) {
      ATCoordinate[0] = (866 * HRAfactor) * iRow; //866m for 1km-HRA; 260m for 300m-HRA
      ATCoordinate[1] = (500 * HRAfactor) * iCol; //500m for 1km-HRA; 150m for 300m-HRA
      ATCoordinate[2] = ICETHICK - ATDepth;
      MirrorATCoordinate[0] = ATCoordinate[0];
      MirrorATCoordinate[1] = ATCoordinate[1];
      MirrorATCoordinate[2] = -(ICETHICK - ATDepth);
   }

   else {
      //get AT's coordinates from the bin row and column number
      ATCoordinate[0] = ATGap * iRow;
      ATCoordinate[1] = ATGap * iCol;
      ATCoordinate[2] = ICETHICK - ATDepth;
      MirrorATCoordinate[0] = ATCoordinate[0];
      MirrorATCoordinate[1] = ATCoordinate[1];
      MirrorATCoordinate[2] = -(ICETHICK - ATDepth);
   }

}

void GetATLocation(int iRow, int iCol, double* ATCoordinate)
{
   if (HEXAGONAL) {
      ATCoordinate[0] = (866 * HRAfactor) * iRow; //866m for 1km-HRA; 260m for 300m-HRA
      ATCoordinate[1] = (500 * HRAfactor) * iCol; //500m for 1km-HRA; 150m for 300m-HRA
      ATCoordinate[2] = ICETHICK - ATDepth;
   }

   else {
      //this is for square arrays
      ATCoordinate[0] = ATGap * iRow;
      ATCoordinate[1] = ATGap * iCol;
      ATCoordinate[2] = ICETHICK - ATDepth;
   }

}

void GetMirrorATLocation(int iRow, int iCol, double* MirrorATCoordinate)
{

   if (HEXAGONAL) {
      MirrorATCoordinate[0] = (866 * HRAfactor) * iRow; //866m for 1km-HRA; 260m for 300m-HRA
      MirrorATCoordinate[1] = (500 * HRAfactor) * iCol; //500m for 1km-HRA; 150m for 300m-HRA
      MirrorATCoordinate[2] = -(ICETHICK - ATDepth);
   }

   else {
//this is for square arrays
      MirrorATCoordinate[0] = ATGap * iRow;
      MirrorATCoordinate[1] = ATGap * iCol;
      MirrorATCoordinate[2] = -(ICETHICK - ATDepth);
   }

}

void GetInteractionPoint(double* posnu)
{

   if (HEXAGONAL) {
      posnu[0] = Rand3.Rndm() * ((2) * 1000. + 2 * MAX_DISTANCE) - MAX_DISTANCE;
      posnu[1] = Rand3.Rndm() * ((2) * 1000. + 2 * MAX_DISTANCE) - MAX_DISTANCE;
      posnu[2] = Rand3.Rndm() * ICETHICK; //postive
   }

//   if(HEXAGONAL){
//    posnu[0]=Rand3.Rndm()*((NROWS-1)*866.+2*EDGE)-EDGE;
//    posnu[1]=Rand3.Rndm()*((NCOLS-1)*500.+2*EDGE)-EDGE;
//    posnu[2]=Rand3.Rndm()*ICETHICK;//postive
//   }


   else {

      //  posnu[0]=Random(0.,101*ATGap );
      // posnu[0]=Random(-4000.,NROWS*ATGap+4000. );
      posnu[0] = Rand3.Rndm() * ((NROWS - 1) * ATGap + 2 * MAX_DISTANCE) - MAX_DISTANCE;
      // cout<<posnu[0]/101./ATGap<<endl;
      // posnu[1]=Random(-4000.,NCOLS*ATGap+4000. );
      posnu[1] = Rand3.Rndm() * ((NCOLS - 1) * ATGap + 2 * MAX_DISTANCE) - MAX_DISTANCE;
      // cout<<posnu[1]/101./ATGap<<endl;
      posnu[2] = Rand3.Rndm() * ICETHICK; //postive
      // cout<<posnu[2]/ICETHICK<<endl;
      // MirrorPoint[0]=posnu[0];
      //MirrorPoint[1]=posnu[1];
      //MirrorPoint[2]=2*ICETHICK-posnu[2];

      // test<<posnu[0]<< "  "<<posnu[1]<<endl;
   }

   //for checking the random number
   // cout<<"rnd1="<<rnd1<<"   rnd2="<<rnd2<who<"  rnd3="<<rnd3<<endl;

}

void GetPosnuBin(double* posnu, int& posnu_iRow, int& posnu_iCol)
{

   if (HEXAGONAL)  {
      if (posnu[0] < 0)
         posnu_iRow = 0;
      else if (posnu[0] >= (NROWS - 1) * (866 * HRAfactor)) //866m for 1km separations; 260m for 300m-HRA
         posnu_iRow = NROWS;
      else
         posnu_iRow = (int)(posnu[0] / (866 * HRAfactor)) + 1;

      if (posnu[1] < 0)
         posnu_iCol = 0;
      else if (posnu[1] >= (NCOLS - 1) * (500 * HRAfactor)) //500m for 1km-HRA; 150m for 300m-HRA
         posnu_iCol = NCOLS;
      else
         posnu_iCol = (int)(posnu[1] / (500 * HRAfactor)) + 1;
   }


   else {
      if (posnu[0] < 0)
         posnu_iRow = 0;
      else if (posnu[0] >= (NROWS - 1)*ATGap)
         posnu_iRow = NROWS;
      else
         posnu_iRow = (int)(posnu[0] / ATGap) + 1;

      if (posnu[1] < 0)
         posnu_iCol = 0;
      else if (posnu[1] >= (NCOLS - 1)*ATGap)
         posnu_iCol = NCOLS;
      else
         posnu_iCol = (int)(posnu[1] / ATGap) + 1;
   }

}

double AttenLengthAtDepth(double d)
{
  d *= 420. / ICETHICK;
  double L = (1250.*0.08886 * exp(-0.048827 * (225.6746 - 86.517596 * log10(848.870 - (d)))));
  L *= ATTEN_UP/262.0;
  return L;
}

void GetAttenlength(double* posnu, double& attenlength_up, double& attenlength_down)
{
   if (CONST_ATTENLENGTH) {
      attenlength_up = ATTEN_UP;
      attenlength_down = ATTEN_DOWN;
      // attenlength_up=400;//m
      // attenlength_down=250.;
   } else {
      double interactiondepth = ICETHICK - posnu[2]; //positive
      double effectivedepth = 0.;
      double effectivemaxdepth = 0.;
      double sum_down = 0;
      double sum_up = 0;
      double xtot = 0;

      int bin;
      double dx = 1; //using 1m as a step
      double dx2 = dx / 2.;

      effectivedepth = interactiondepth;
      effectivemaxdepth = ICETHICK;

      //cout << "ICETHICK = "<<ICETHICK<<endl;
      //cout << "z = "<<posnu[2]<<endl;
      //cout<<"interactiondepth = "<<interactiondepth<<endl;
      //cout << "effectivedepth = "<<effectivedepth<<endl;

      dx = (effectivemaxdepth - effectivedepth) * 0.005;
      if (dx < 1.e-5)
         dx = 1.e-5;
      dx2 = 0.5 * dx;

      for (bin = 0; bin < (round)((effectivemaxdepth - effectivedepth) / dx); bin++){
        	sum_down += dx * AttenLengthAtDepth(dx * bin + dx2 + effectivedepth); //down to the ice bottom
        	xtot += dx;
	       }


      dx = effectivemaxdepth * 0.005;
      if (dx < 1.e-5)
         dx = 1.e-5;
      dx2 = 0.5 * dx;

      for (bin = 0; bin < (round)((effectivemaxdepth) / dx); bin++){
	       sum_down += dx * AttenLengthAtDepth(dx * bin + dx2); //upward to the surface
	       xtot += dx;
        }

      //cout << "upward: xtot = " << xtot << endl;
      //cout << "upward: 2*effectivemaxdepth - effectivedepth = " << (2 * effectivemaxdepth - effectivedepth) << endl<<endl;

      if (effectivedepth < 1.e-1)
         effectivedepth = 1.e-1;

      dx = effectivedepth * 0.005;
      if (dx < 1.e-5)
         dx = 1.e-5;
      dx2 = dx / 2.;

      for (bin = 0; bin < (round)(effectivedepth / dx); bin++) {
	       sum_up += dx * AttenLengthAtDepth(dx * bin + dx2);
      }
      attenlength_up = sum_up / effectivedepth;
      attenlength_up *= ATTEN_FACTOR;

   }

   //lamda_up<<ICETHICK-posnu[2]<<"  "<<attenlength_up<<endl;
   //lamda_down<<ICETHICK-posnu[2]<<"  "<<attenlength_down<<endl;

   //  cout<<"attenup="<<attenlength_up<<" attendown="<<attenlength_down<<endl;
   //  test<<posnu[2]<<"   "<<attenlength<<endl;



}


void GetEntryPoint(double theta_nu, double* nnu, double* posnu, double* entry) //get ice entry point
{
   if (theta_nu > PI / 2 && theta_nu <= PI) { //downward
      entry[0] = -fabs((ICETHICK - posnu[2]) / nnu[2]) * nnu[0] + posnu[0];
      entry[1] = -fabs((ICETHICK - posnu[2]) / nnu[2]) * nnu[1] + posnu[1];
      entry[2] = -fabs((ICETHICK - posnu[2]) / nnu[2]) * nnu[2] + posnu[2];
   } else if (theta_nu <= PI / 2) { //upward
      entry[0] = posnu[0] - nnu[0] * posnu[2] / nnu[2];
      entry[1] = posnu[1] - nnu[1] * posnu[2] / nnu[2];
      entry[2] = posnu[2] - nnu[2] * posnu[2] / nnu[2]; //nnu[2]=cos(theta_nu);
      // cout<<theta_nu<<"  =?"<<acos(nnu[2])<<endl;
   } else
      cout << "GetEntryPoint: wrong theta_nu [" << theta_nu << "]" << endl;

}


void GetBouncePoint(double theta_nposnu2MirrorAT, double* nposnu2MirrorAT, double* posnu, double* bounce)
//KD:06/29/11, get water bounce point for reflected
{

   if (theta_nposnu2MirrorAT > PI / 2 && theta_nposnu2MirrorAT <= PI) { //downward
      bounce[0] = fabs((0 - posnu[2]) / nposnu2MirrorAT[2]) * nposnu2MirrorAT[0] + posnu[0];
      bounce[1] = fabs((0 - posnu[2]) / nposnu2MirrorAT[2]) * nposnu2MirrorAT[1] + posnu[1];
      bounce[2] = fabs((0 - posnu[2]) / nposnu2MirrorAT[2]) * nposnu2MirrorAT[2] + posnu[2];
   }

   /*
     else if(theta_nu<=PI/2)//upward
       {
         entry[0]=posnu[0]-nnu[0]*posnu[2]/nnu[2];
         entry[1]=posnu[1]-nnu[1]*posnu[2]/nnu[2];
         entry[2]=posnu[2]-nnu[2]*posnu[2]/nnu[2];//nnu[2]=cos(theta_nu);
         // cout<<theta_nu<<"  =?"<<acos(nnu[2])<<endl;
       }
     */
   /*
   else
     cout<<"GetBouncePoint: wrong theta_nposnu2MirrorAT ["<<theta_nposnu2MirrorAT
         << "] bounce=("<<bounce[0]<<", "<<bounce[1]<<", "<<bounce[2]<")"<<endl;
   */

}


/*double  Random(double n1, double n2)
  {
  if(n1>n2) return -1;
  if(n1==n2) return 0;
  return n1+double((n2-n1))*double(rand())/RAND_MAX;

  }*/

double Distance(double* point1, double* point2)
{
   return sqrt(Square(point1[0] - point2[0]) + Square(point1[1] - point2[1]) + Square(point1[2] - point2[2]));



}

double Square(double a)
{
   return a * a;

}

double GetDecayLength(double EXPONENT)
{
   //Gives Tau decay lengeth in meters
   //Cut off length at 10^18 eV, since Tau energy loss is not yet modeled.
   //This will reduce the expected decay length at high energies.
   if (EXPONENT>18) {
	return 50.0*pow(10,3);
   }
   else {
      return 50.0*pow(10,EXPONENT-15);
   }
}

void Sigma(double pnu, double& sigma, double& int_len_kgm2)
{
   sigma = (7.84E-40) * pow(pnu / 1.E9, 0.363) * SIGMA_FACTOR; //m^2
   int_len_kgm2 = AMU / sigma;
   //  return sigma;


}


/*double VmMHz1m(double freq,double pnu)//get the R|E| in changle direction
//R*|(vector)E(w,R,changle)|
//the unit is V/MHz,electric field(V/m) multiplied by observation distance(1m)
{
double vmmhz1m;
double freq0=1150.;//MHz
return vmmhz1m=2.53E-7*(pnu/1E12)*(freq/freq0)*(1/(1+pow(freq/freq0,1.44)));
//astro-ph/0003315 v1 22 Mar 2000

}*/

void GetEmHadFrac(string nuflavor, string current, double elast_y, double theta_nu, double L_TauDecay, double L_Ice, double& emfrac, double& hadfrac)
{

   if (nuflavor == "nue" && current == "cc") {
      emfrac = 1. - elast_y;
      hadfrac = elast_y;
   } else if (nuflavor == "numu" && current == "cc") {
      emfrac = 1.E-10;
      hadfrac = elast_y;
   } else if (nuflavor == "nutau" && current == "cc") {
      // FW: behaves like a muon
      // KD: Now changed to SB
      emfrac = 1.E-10;
      //down going, with double bang in ice
      if (theta_nu >= PI/2.0 && L_Ice > L_TauDecay) {
         if (elast_y < 0.5)
            hadfrac = 1. - elast_y;
         else
            hadfrac = elast_y;
      }
      //down going with double bang below ice
      else if (theta_nu >= PI/2.0 && L_Ice <= L_TauDecay) {
            hadfrac = elast_y;
      }
      //up going
      else if (theta_nu < PI/2.0) {
         if (elast_y < 0.5)
            hadfrac = 1. - elast_y;
         else
            hadfrac = elast_y;
      }
   }


   else if (current == "nc") {
      emfrac = 1.E-10;
      hadfrac = elast_y;
   }
   //GetCurrent(current);
   /* if(nuflavor=="nue")
      {emfrac=0.8;
      hadfrac=0.2;
      }
      else if(nuflavor=="numu"||nuflavor=="nutau")
      {emfrac=0.;
      hadfrac=0.2;
      }*/
   // emfrac=1.;
   // hadfrac=0.;

}

double VmMHz1m(double viewangle, double freq, double emfrac, double hadfrac, double em_shower_length, double had_shower_length) //using the pulse scaling of David Seckel
{
   double nuy = 0; //MHz
   double em_nuz = 0;
   double had_nuz = 0;
   double Gy = 0;
   double em_Gz = 0;
   double had_Gz = 0;

   double delta = fabs(viewangle - changle); //steradian


   // double f0=2.53E-7/freq0/sin(changle);
   // cout<<"sing(changle)"<<sin(changle)<<endl;

   // cout<<"f0="<<f0<<endl;

   // cout<<had_shower_length<<endl; //so far it's correct
   nuy = freq * (1 + delta / tan(changle));
   //  cout<<"freq="<<freq<<"nuy="<<nuy<<endl;
   em_nuz = freq * em_shower_length * delta;
   had_nuz = freq * had_shower_length * delta;
   // cout<<freq<<"  "<<em_nuz<<endl;

   Gy = 1 / (1 + pow(nuy / freq0, 1.44));
   // cout<<freq<<"  "<<Gy<<endl;
   // em_Gz=exp(-1.*pow(em_nuz/1000/0.032,2));
   em_Gz = exp(-1.*pow(em_nuz / 32., 2));
   //had_Gz=exp(-1.*pow(had_nuz/1000/0.032,2));
   had_Gz = exp(-1.*pow(had_nuz / 32., 2));

   double pnu_12 = pnu / 1.E12;
   // cout<<"em_Gz"<<em_Gz<<endl;
   double em_vmmhz1m = f0 * freq * sin(viewangle) * Gy * em_Gz * pnu_12;
   double had_vmmhz1m = f0 * freq * sin(viewangle) * Gy * had_Gz * pnu_12;
   double vmmhz1m = emfrac * em_vmmhz1m + hadfrac * had_vmmhz1m;
   //  cout<<em_vmmhz1m<<"   "<<had_vmmhz1m<<endl;
   //  cout<<emfrac<< "   "<<hadfrac<<endl;
   //cout<<"f0="<<f0<<endl;
   // cout<<"freq="<<freq<<endl;
   // cout<<"viewangle"<<viewangle*RAD2DEG<<endl;
   //cout<<"pnu"<<pnu<<endl;
   //  cout<<"em_Gy"<<Gy<<endl;

   //  cout<<vmmhz1m<<endl;
   return vmmhz1m;

}

/*double VmMHz(double vmmhz1m, double *posnu,double *ATCoordinate)
  {

  return vmmhz1m/Distance(posnu,ATCoordinate);//electric field


  }*/
double VmMHz(double vmmhz1m, double r)
{
   return vmmhz1m / r; //electric field
}


void Times(double a, double* v1, double* v2)
{
   for (int i = 0; i < 3; i++) {
      v2[i] = a * v1[i];

   }


}

double Dot(double* v1, double* v2)
{
   return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
void Cross(double* v1, double* v2, double* v3)
{
   v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
   v3[1] = -1 * v1[0] * v2[2] + v1[2] * v2[0];
   v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

double Mag(double* v1)//get the magnitude of a vector
{
   return sqrt(Square(v1[0]) + Square(v1[1]) + Square(v1[2]));
}

void  nVector(double* v1, double* nv1)
{
   nv1[0] = v1[0] / Mag(v1);
   nv1[1] = v1[1] / Mag(v1);
   nv1[2] = v1[2] / Mag(v1);
}

void outVector(double* v, int n)
{
   for (int i = 0; i < n; i++) {
      cout << v[i] << "   ";
   }
   cout << endl;

}

//KD to extract spherical angles from a given normalized vector
// Yeah, yeah, I can be more sophisticated and use Vector's all over the place, but no time to implement that
void SphAngles(double* v1, double* sphangles)
{
//double sphangles[2];
   if (v1[0] > 0) {
      if (v1[1] > 0)
         sphangles[0] = atan(v1[1] / v1[0]);
      if (v1[1] < 0)
         sphangles[0] = 2 * PI + atan(v1[1] / v1[0]);
   } else if (v1[0] < 0) {
      if (v1[1] > 0) {
         sphangles[0] = PI + atan(v1[1] / v1[0]);
      }
      if (v1[1] < 0) {
         sphangles[0] = PI + atan(v1[1] / v1[0]);
      }
   }
   sphangles[1] = acos(v1[2]);
//I have to improve upon this above!
}


double Angle(double* v1, double* v2)
{
   return acos(Dot(v1, v2) / Mag(v1) / Mag(v2));
}

void VectorMinus(double* v1, double* v2, double* v3)
{
   v3[0] = v1[0] - v2[0];
   v3[1] = v1[1] - v2[1];
   v3[2] = v1[2] - v2[2];
}

void VectorPlus(double* v1, double* v2, double* v3)
{
   for (int i = 0; i < 3; i++) {
      v3[i] = v1[i] + v2[i];
   }
}

double Min(double* array, int n)
{
   double min = 1.E24;
   for (int k = 0; k < n; k++) {
      if (array[k] < min)
         min = array[k];
   }
   return min;
}



double Max(double* array, int n)
{
   double max = -1.E24;
   for (int k = 0; k < n; k++) {
      if (array[k] > max)
         max = array[k];
   }
   return max;

}


/*int GetEm()
  {
  int em;
  double temp=Random(0,1);
  if(temp<0.63212056)
  em=1;
  else
  em=0;
  return em;

  }*/

double VmMHz_attenuated(double d_posnu2AT, double vmmhz, double attenlength) //needing modification
{
   double temp = d_posnu2AT / attenlength;
   if (temp < 20) {
      return vmmhz * exp(-temp);
   } else
      return 0;
}


double GetNoise(double bw)
{
   //the unit of bw is MHz
   //double TICE = 245.; //K ice noise
   //double TSYS = 200.; //K system noise

   double TNOISE = 350;
   double vnoise = sqrt(TNOISE * 50.*KBOLTZ * bw * 1.E6);
   return vnoise;
}

void GetChord(double theta_nu, double& chord_inCRUST, double& chord2_inICE)
{
   double tangle = PI / 2. - theta_nu; //angle between nnu and CRUST surface
   chord_inCRUST = 2 * R_earth * sin(tangle);
   chord2_inICE = ICETHICK / sin(tangle);
}


double GetWeight(double chord_inICE, double chord_inCRUST, double chord2_inICE, double sigma, double theta_nu)
{
   double temp1 = exp(-chord_inICE * sigma * DensityICE / AMU);
   double temp2 = exp(-chord_inCRUST * sigma * DensityCRUST / AMU);
   double temp3 = exp(-chord2_inICE * sigma * DensityICE / AMU);
   if (theta_nu > PI / 2)
      return temp1;
   else
      return temp1 * temp2 * temp3;
}


double GetTauRegen(string current, double energy, double theta_nu, double L_TauDecay, double H)
// KD 10/2010
{
   double tauweight = 0;
   int NEbins = 0;
// data below for GZK spectrum
//if (GZK){ JET:  If TauRegeneration is on, we want GZK table to be used
// regardless of wether or not GZK is toggled on.
   double tautable[10][9];
   double angle_range[10];
   double energy_range[11];

   double theta_nu_deg = theta_nu * RAD2DEG;
   double D = ICETHICK - H;
   double GeoFactor = 1.0;
   //Add to effective volume for upward going Tau, since original interaction can occur below the ice
   if (current == "cc" && (L_TauDecay*cos(theta_nu) + D) > ICETHICK){
      GeoFactor = (L_TauDecay*cos(theta_nu) + D)/ICETHICK;
   }

   if (GZK){ //GZK Spectrum
     NEbins=10;
  // energy 15.0-15.5
     tautable[0][0] = 3.58 ;
     tautable[0][1] = 3.01 ;
     tautable[0][2] = 2.38 ;
     tautable[0][3] = 1.94 ;
     tautable[0][4] = 1.63 ;
     tautable[0][5] = 1.43 ;
     tautable[0][6] = 1.29 ;
     tautable[0][7] = 1.16 ;
     tautable[0][8] = 1.05 ;

  // energy 15.5-16.0
     tautable[1][0] = 2.92;
     tautable[1][1] = 2.94;
     tautable[1][2] = 2.65;
     tautable[1][3] = 2.28;
     tautable[1][4] = 1.87;
     tautable[1][5] = 1.74 ;
     tautable[1][6] = 1.47;
     tautable[1][7] = 1.27;
     tautable[1][8] = 1.07;

  // energy 16.0-16.5
     tautable[2][0] = 1.33;
     tautable[2][1] = 1.89;
     tautable[2][2] = 2.23;
     tautable[2][3] = 2.22;
     tautable[2][4] = 2.11;
     tautable[2][5] = 1.89;
     tautable[2][6] = 1.69;
     tautable[2][7] = 1.42;
     tautable[2][8] = 1.13;

  // energy 16.5-17.0
     tautable[3][0] = 0.26;
     tautable[3][1] = 0.57;
     tautable[3][2] = 0.97;
     tautable[3][3] = 1.27;
     tautable[3][4] = 1.46;
     tautable[3][5] = 1.51;
     tautable[3][6] = 1.48;
     tautable[3][7] = 1.37;
     tautable[3][8] = 1.14;

  // energy 17.0-17.5
     tautable[4][0] = 0.02;
     tautable[4][1] = 0.06;
     tautable[4][2] = 0.17;
     tautable[4][3] = 0.34;
     tautable[4][4] = 0.52;
     tautable[4][5] = 0.65;
     tautable[4][6] = 0.80;
     tautable[4][7] = 0.92;
     tautable[4][8] = 0.98;

  // energy 17.5-18.0
     tautable[5][0] = 0.00;
     tautable[5][1] = 0.00;
     tautable[5][2] = 0.02;
     tautable[5][3] = 0.08;
     tautable[5][4] = 0.17;
     tautable[5][5] = 0.28;
     tautable[5][6] = 0.43;
     tautable[5][7] = 0.64;
     tautable[5][8] = 0.87;

  // energy 18.0-18.5
     tautable[6][0] = 0.00;
     tautable[6][1] = 0.00;
     tautable[6][2] = 0.00;
     tautable[6][3] = 0.01;
     tautable[6][4] = 0.03;
     tautable[6][5] = 0.07;
     tautable[6][6] = 0.17;
     tautable[6][7] = 0.36;
     tautable[6][8] = 0.77;

  // energy 18.5-19.0
     tautable[7][0] = 0.00;
     tautable[7][1] = 0.00;
     tautable[7][2] = 0.00;
     tautable[7][3] = 0.00;
     tautable[7][4] = 0.00;
     tautable[7][5] = 0.01;
     tautable[7][6] = 0.05;
     tautable[7][7] = 0.22;
     tautable[7][8] = 0.74;

  // energy 19.0-19.5
     tautable[8][0] = 0.00;
     tautable[8][1] = 0.00;
     tautable[8][2] = 0.00;
     tautable[8][3] = 0.00;
     tautable[8][4] = 0.00;
     tautable[8][5] = 0.00;
     tautable[8][6] = 0.01;
     tautable[8][7] = 0.15;
     tautable[8][8] = 0.73;

  // energy 19.5-20.0
     tautable[9][0] = 0.00;
     tautable[9][1] = 0.00;
     tautable[9][2] = 0.00;
     tautable[9][3] = 0.00;
     tautable[9][4] = 0.00;
     tautable[9][5] = 0.00;
     tautable[9][6] = 0.02;
     tautable[9][7] = 0.13;
     tautable[9][8] = 0.67;

     //angle Range
     angle_range[0] = 60.;
     angle_range[1] = 65.;
     angle_range[2] = 70.;
     angle_range[3] = 75.;
     angle_range[4] = 80.;
     angle_range[5] = 82.;
     angle_range[6] = 84.;
     angle_range[7] = 86.;
     angle_range[8] = 88.;
     angle_range[9] = 90.;

     //energy_range
     energy_range[0] = 15.0 ;
     energy_range[1] = 15.5 ;
     energy_range[2] = 16.0 ;
     energy_range[3] = 16.5 ;
     energy_range[4] = 17.0 ;
     energy_range[5] = 17.5 ;
     energy_range[6] = 18.0 ;
     energy_range[7] = 18.5 ;
     energy_range[8] = 19.0 ;
     energy_range[9] = 19.5 ;
     energy_range[10] = 20.0;
  }// end of if GZK
  else{ //E^2 Spectrum
     NEbins=6;
    // energy 15.0-16.0
     tautable[0][0]=0.59 ;
     tautable[0][1]=0.69 ;
     tautable[0][2]=0.78 ;
     tautable[0][3]=0.84 ;
     tautable[0][4]=0.89 ;
     tautable[0][5]=0.90 ;
     tautable[0][6]=0.93 ;
     tautable[0][7]=0.96 ;
     tautable[0][8]=0.99 ;

    // energy 16.0-17.0
     tautable[1][0]=0.19;
     tautable[1][1]=0.33;
     tautable[1][2]=0.49;
     tautable[1][3]=0.60;
     tautable[1][4]=0.70;
     tautable[1][5]=0.79;
     tautable[1][6]=0.84;
     tautable[1][7]=0.91;
     tautable[1][8]=0.97;

    // energy 17.0-18.0
     tautable[2][0]=0.01;
     tautable[2][1]=0.03;
     tautable[2][2]=0.08;
     tautable[2][3]=0.22;
     tautable[2][4]=0.35;
     tautable[2][5]=0.39;
     tautable[2][6]=0.66;
     tautable[2][7]=0.75;
     tautable[2][8]=0.93;

    // energy 18.0-19.0
     tautable[3][0]=0.00;
     tautable[3][1]=0.00;
     tautable[3][2]=0.00;
     tautable[3][3]=0.01;
     tautable[3][4]=0.04;
     tautable[3][5]=0.07;
     tautable[3][6]=0.15;
     tautable[3][7]=0.36;
     tautable[3][8]=0.77;

    // energy 19.0-20.0
     tautable[4][0]=0.00;
     tautable[4][1]=0.00;
     tautable[4][2]=0.00;
     tautable[4][3]=0.00;
     tautable[4][4]=0.00;
     tautable[4][5]=0.00;
     tautable[4][6]=0.01;
     tautable[4][7]=0.21;
     tautable[4][8]=0.73;

    // energy 20.0-21.5
     tautable[5][0]=0.00;
     tautable[5][1]=0.00;
     tautable[5][2]=0.00;
     tautable[5][3]=0.00;
     tautable[5][4]=0.00;
     tautable[5][5]=0.00;
     tautable[5][6]=0.00;
     tautable[5][7]=0.35;
     tautable[5][8]=0.90;

     //angle range
     angle_range[0]=60.;
     angle_range[1]=65.;
     angle_range[2]=70.;
     angle_range[3]=75.;
     angle_range[4]=80.;
     angle_range[5]=82.;
     angle_range[6]=84.;
     angle_range[7]=86.;
     angle_range[8]=88.;
     angle_range[9]=90.;

     //energy range
     energy_range[0]=15.0 ;
     energy_range[1]=16.0 ;
     energy_range[2]=17.0 ;
     energy_range[3]=18.0 ;
     energy_range[4]=19.0 ;
     energy_range[5]=20.0 ;
     energy_range[6]=21.5 ;

  }//end else

//double tauweight=0;

   for (int j = 0; j < NEbins; j++) {
      if (log10(energy) >= energy_range[j] && log10(energy) < energy_range[j + 1]) {
         for (int i = 0; i < 9; i++) {
            if (theta_nu_deg >= angle_range[i] && theta_nu_deg < angle_range[i + 1]) {
               tauweight = tautable[j][i] * GeoFactor;
            }
         }
      }
   }

   return tauweight;
}


void GetMaxDis(double pnu, double& Max_distance)
{
   if (pnu < 1e19)
      Max_distance = 3000.;
   else if (pnu < 1e21)
      Max_distance = 4000;
   else
      Max_distance = 10000; //meters



}

TStyle* RootStyle()
{


   TStyle* RootStyle = new TStyle("Root-Style", "The Perfect Style for Plots ;-)");

#ifdef __CINT__
   TStyle* GloStyle;
   GloStyle = gStyle;                          // save the global style reference

   gStyle = RootStyle;
#endif
   // otherwise you need to call TROOT::SetStyle("Root-Style")

   // Paper size

   RootStyle->SetPaperSize(TStyle::kUSLetter);

   // Canvas

   RootStyle->SetCanvasColor(0);
   RootStyle->SetCanvasBorderSize(10);
   RootStyle->SetCanvasBorderMode(0);
   RootStyle->SetCanvasDefH(600);
   RootStyle->SetCanvasDefW(600);
   RootStyle->SetCanvasDefX(10);
   RootStyle->SetCanvasDefY(10);

   // Pads

   RootStyle->SetPadColor(0);
   RootStyle->SetPadBorderSize(10);
   RootStyle->SetPadBorderMode(0);
   //  RootStyle->SetPadBottomMargin(0.13);
   RootStyle->SetPadBottomMargin(0.16);
   RootStyle->SetPadTopMargin(0.08);
   RootStyle->SetPadLeftMargin(0.18);
   RootStyle->SetPadRightMargin(0.05);
   RootStyle->SetPadGridX(0);
   RootStyle->SetPadGridY(0);
   RootStyle->SetPadTickX(1);
   RootStyle->SetPadTickY(1);

   // Frames

   RootStyle->SetFrameFillStyle(0);
   RootStyle->SetFrameFillColor(0);
   RootStyle->SetFrameLineColor(1);
   RootStyle->SetFrameLineStyle(0);
   RootStyle->SetFrameLineWidth(2);
   RootStyle->SetFrameBorderSize(10);
   RootStyle->SetFrameBorderMode(0);


   // Histograms

   RootStyle->SetHistFillColor(0);
   RootStyle->SetHistFillStyle(1);
   RootStyle->SetHistLineColor(1);
   RootStyle->SetHistLineStyle(0);
   RootStyle->SetHistLineWidth(3);

   // Functions

   RootStyle->SetFuncColor(1);
   RootStyle->SetFuncStyle(0);
   RootStyle->SetFuncWidth(1);

   //Legends

   RootStyle->SetStatBorderSize(2);
   RootStyle->SetStatFont(42);
   //  RootStyle->SetOptStat       (111111);
   RootStyle->SetOptStat(0);
   RootStyle->SetStatColor(0);
   RootStyle->SetStatX(0.93);
   RootStyle->SetStatY(0.90);
   RootStyle->SetStatFontSize(0.07);
   //  RootStyle->SetStatW         (0.2);
   //  RootStyle->SetStatH         (0.15);

   // Labels, Ticks, and Titles

   RootStyle->SetTickLength(0.015, "X");
   RootStyle->SetTitleSize(0.10, "X");
   RootStyle->SetTitleOffset(0.800, "X");
   RootStyle->SetTitleBorderSize(0);
   //  RootStyle->SetTitleFontSize((double)3.);
   RootStyle->SetLabelOffset(0.015, "X");
   RootStyle->SetLabelSize(0.050, "X");
   RootStyle->SetLabelFont(42   , "X");

   RootStyle->SetTickLength(0.015, "Y");
   RootStyle->SetTitleSize(0.10, "Y");
   RootStyle->SetTitleOffset(0.600, "Y");
   RootStyle->SetLabelOffset(0.015, "Y");
   RootStyle->SetLabelSize(0.050, "Y");
   RootStyle->SetLabelFont(42   , "Y");

   RootStyle->SetTitleFont(42, "XY");
   RootStyle->SetTitleColor(1);

   // Options

   RootStyle->SetOptFit(1);

   RootStyle->SetMarkerStyle(20);
   RootStyle->SetMarkerSize(0.4);

   //  cout << ">> Style initialized with the Root Style!" << endl;
   //  cout << ">> " << modified << endl << endl;
   return RootStyle;
}

void CloseTFile(TFile* hfile)
{
   hfile->Write();
   hfile->Close();
}

double GetLPM(double X0DEPTH)
{
   return 2.E15 * (X0DEPTH / X0ICE); //2Pev LPM enerny
}


double Gety()
{
   float rnd;
   float x = 0;
   const double R1 = 0.36787944;
   const double R2 = 0.63212056;
   rnd = Rand3.Rndm(1); //(0,1)
   x = TMath::Power(-log(R1 + rnd * R2), 2.5);
   return x;
}

/*
int GetBeamWidths(double flare[4][NFREQ],double gain[2][NFREQ],double freq[NFREQ]){
//double GetGainquick(double *flare_v, double *flare_h, double* gainv){
//flare_v = (30416/gainv)*DEG2RAD;
//flare_h = (114927/gainv)*DEG2RAD;

  return 1;
}
*/

void GetFlare(double freq, double* flare_tmp)
{
   double specs[5][4];
   double freq_specs[5];
   double scale = 0;

//KD: here are the beam widths from Jordan, but flip it due to hitangle e and h definitions in main code
   //GOOD ONES
   specs[0][0] = 144;
   specs[0][1] = 74;
   specs[0][2] = 0;
   specs[0][3] = 0;
   specs[1][0] = 128;
   specs[1][1] = 70;
   specs[1][2] = 0;
   specs[1][3] = 0;
   specs[2][0] = 134;
   specs[2][1] = 74;
   specs[2][2] = 0;
   specs[2][3] = 0;
   specs[3][0] = 132;
   specs[3][1] = 68;
   specs[3][2] = 0;
   specs[3][3] = 0;
   specs[4][0] = 154;
   specs[4][1] = 70;
   specs[4][2] = 0;
   specs[4][3] = 0;

   freq_specs[0] = 300; //  freq_specs[0]=300.E6;
   freq_specs[1] = 600; //  freq_specs[1]=600.E6;
   freq_specs[2] = 900; //  freq_specs[2]=900.E6;
   freq_specs[3] = 1200; //  freq_specs[3]=1200.E6;
   freq_specs[4] = 1500; //  freq_specs[4]=1500.E6;

   if (freq < freq_specs[0]) {
     for (int j = 0; j < 4; j++) {
       flare_tmp[j] = specs[0][j] * DEG2RAD;
     } //for
   } //if
   else if (freq >= freq_specs[3]) {
     for (int j = 0; j < 4; j++) {
       flare_tmp[j] = specs[3][j] * DEG2RAD;
     } //for
   } //else if
   else {
     for (int i = 0; i < 4; i++) {
       if (freq >= freq_specs[i] && freq < freq_specs[i + 1]) {
	 scale = (freq - freq_specs[i]) / (freq_specs[i + 1] - freq_specs[i]);

	 for (int j = 0; j < 4; j++) {
	   flare_tmp[j] = (specs[i][j] + scale * (specs[i + 1][j] - specs[i][j])) * DEG2RAD;
	 } //for
	 i = 4;
       } //if
     } //for
   } //else

}

int GetBeamWidths(double flare[4][NFREQ], double gain[2][NFREQ], double freq[NFREQ])
{

   // first component is frequency
   // second component is which plane and which polarization
   // it goes  e-plane: vp/hp, h-plane: vp/hp

   // these number were read from antenna specs

   double specs[5][4];
   /*
   //KD: below is the original specs for the horn
     specs[0][0]=57.5;  specs[0][1]=58.5;  specs[0][2]=66;  specs[0][3]=57;
     specs[1][0]=33.5;  specs[1][1]=34.5;  specs[1][2]=36.5;  specs[1][3]=38;
     specs[2][0]=50.5;  specs[2][1]=53;  specs[2][2]=33;  specs[2][3]=32;
     specs[3][0]=43.5;  specs[3][1]=43;  specs[3][2]=39;  specs[3][3]=41.5;
     specs[4][0]=36.5;  specs[4][1]=46.5;  specs[4][2]=32;  specs[4][3]=31;
   */


//KD: here are the beam widths from Jordan, but flip it due to hitangle e and h definitions in main code
   //GOOD ONES
   specs[0][0] = 144;
   specs[0][1] = 74;
   specs[0][2] = 0;
   specs[0][3] = 0;
   specs[1][0] = 128;
   specs[1][1] = 70;
   specs[1][2] = 0;
   specs[1][3] = 0;
   specs[2][0] = 134;
   specs[2][1] = 74;
   specs[2][2] = 0;
   specs[2][3] = 0;
   specs[3][0] = 132;
   specs[3][1] = 68;
   specs[3][2] = 0;
   specs[3][3] = 0;
   specs[4][0] = 154;
   specs[4][1] = 70;
   specs[4][2] = 0;
   specs[4][3] = 0;

   /* OLD ONES
     specs[0][0]=74;  specs[0][1]=144;  specs[0][2]=0;  specs[0][3]=0;
     specs[1][0]=70;  specs[1][1]=128;  specs[1][2]=0;  specs[1][3]=0;
     specs[2][0]=74;  specs[2][1]=134;  specs[2][2]=0;  specs[2][3]=0;
     specs[3][0]=68;  specs[3][1]=132;  specs[3][2]=0;  specs[3][3]=0;
     specs[4][0]=70;  specs[4][1]=154;  specs[4][2]=0;  specs[4][3]=0;
   */


   /*
   //KD: averaging above as per SB , but that's not what we meant I think
     specs[0][0]=109;  specs[0][1]=109;  specs[0][2]=0;  specs[0][3]=0;
     specs[1][0]=99;  specs[1][1]=99;  specs[1][2]=0;  specs[1][3]=0;
     specs[2][0]=104;  specs[2][1]=104;  specs[2][2]=0;  specs[2][3]=0;
     specs[3][0]=100;  specs[3][1]=100;  specs[3][2]=0;  specs[3][3]=0;
     specs[4][0]=112;  specs[4][1]=112;  specs[4][2]=0;  specs[4][3]=0;
   */

//KD: here are quick beam widths from relationship
// numbers based on taking gain=8 as basis
// specs[0][0]=sqrt(43808/gainv);  specs[0][1]=sqrt(165888/gainv);  specs[0][2]=0;  specs[0][3]=0;
   /*
   specs[0][0]=sqrt(43808/gainv);  specs[0][1]=sqrt(165888/gainv);  specs[0][2]=0;  specs[0][3]=0;

   specs[1][0]=0;  specs[1][1]=0;  specs[1][2]=0;  specs[1][3]=0;
   specs[2][0]=0;  specs[2][1]=0;  specs[2][2]=0;  specs[2][3]=0;
   specs[3][0]=0;  specs[3][1]=0;  specs[3][2]=0;  specs[3][3]=0;
   specs[4][0]=0;  specs[4][1]=0;  specs[4][2]=0;  specs[4][3]=0;

   specs[1][0]=sqrt(43808/gainv);  specs[1][1]=128;  specs[1][2]=0;  specs[1][3]=0;
   specs[2][0]=sqrt(43808/gainv);  specs[2][1]=134;  specs[2][2]=0;  specs[2][3]=0;
   specs[3][0]=sqrt(43808/gainv);  specs[3][1]=132;  specs[3][2]=0;  specs[3][3]=0;
   specs[4][0]=sqrt(43808/gainv);  specs[4][1]=154;  specs[4][2]=0;  specs[4][3]=0;
   */

   /*specs[0][0]=sqrt(41253/gainv); specs[0][1]=sqrt(41253/gainv); specs[0][2]=0;  specs[0][3]=0;
   specs[1][0]=sqrt(41253/gainv); specs[1][1]=sqrt(41253/gainv); specs[1][2]=0;  specs[1][3]=0;
   specs[2][0]=sqrt(41253/gainv); specs[2][1]=sqrt(41253/gainv); specs[2][2]=0;  specs[2][3]=0;
   specs[3][0]=sqrt(41253/gainv); specs[3][1]=sqrt(41253/gainv); specs[3][2]=0;  specs[3][3]=0;
   specs[4][0]=sqrt(41253/gainv); specs[4][1]=sqrt(41253/gainv); specs[4][2]=0;  specs[4][3]=0;
   */

   double specs2[5][2];
   //KD: original Gains specs for horn antennas?
//  specs2[0][0]=8.5;  specs2[0][1]=8.8;
//  specs2[1][0]=11.0;  specs2[1][1]=9.2;
//  specs2[2][0]=9.3;  specs2[2][1]=9.6;
//  specs2[3][0]=10.1;  specs2[3][1]=11.5;
//  specs2[4][0]=8.9;  specs2[4][1]=9.0;
   //--------------------------------------

//KD: here are the power ratios I derived from Jordan's input, note they are NOT gain in dB!
//KD: the specs2 values go in to gain calculation, but these are not the 'gains' used. We just use the 'flare' from here.
//KD: those gain changes should go hand in hand with flare changes
   specs2[0][0] = 6.0;
   specs2[0][1] = 2.5;
   specs2[1][0] = 6.8;
   specs2[1][1] = 2.6;
   specs2[2][0] = 4.2;
   specs2[2][1] = 2.7;
   specs2[3][0] = 6.2;
   specs2[3][1] = 2.5;
   specs2[4][0] = 6.2;
   specs2[4][1] = 2.5;


//KD: 08/27 It looks like FW had the freq in MHz units, and not Hz, so the loop was stuck on freq_specs[0] for all the range 100-1000.
   double freq_specs[5];
   freq_specs[0] = 300; //  freq_specs[0]=300.E6;
   freq_specs[1] = 600; //  freq_specs[1]=600.E6;
   freq_specs[2] = 900; //  freq_specs[2]=900.E6;
   freq_specs[3] = 1200; //  freq_specs[3]=1200.E6;
   freq_specs[4] = 1500; //  freq_specs[4]=1500.E6;

   double scale = 0;

   for (int k = 0; k < NFREQ; k++) {

      if (freq[k] < freq_specs[0]) {
         for (int j = 0; j < 4; j++) {
            flare[j][k] = specs[0][j] * DEG2RAD;
         } //for
      } //if
      else if (freq[k] >= freq_specs[3]) {
         for (int j = 0; j < 4; j++) {
            flare[j][k] = specs[3][j] * DEG2RAD;
         } //for
      } //else if
      else {
         for (int i = 0; i < 4; i++) {
            if (freq[k] >= freq_specs[i] && freq[k] < freq_specs[i + 1]) {
               scale = (freq[k] - freq_specs[i]) / (freq_specs[i + 1] - freq_specs[i]);

               for (int j = 0; j < 4; j++) {
                  flare[j][k] = (specs[i][j] + scale * (specs[i + 1][j] - specs[i][j])) * DEG2RAD;
               } //for
               i = 4;
            } //if
         } //for
      } //else
   } //for (frequencies)


   for (int k = 0; k < NFREQ; k++) {

      if (freq[k] < freq_specs[0]) {
         for (int j = 0; j < 2; j++) {
            gain[j][k] = 10 * log10(specs2[0][j]);
         } //for
      } //if
      else if (freq[k] >= freq_specs[3]) {
         for (int j = 0; j < 2; j++) {
            gain[j][k] = 10 * log10(specs2[3][j]);
         } //for
      } //else if
      else {
         for (int i = 0; i < 4; i++) {
            if (freq[k] >= freq_specs[i] && freq[k] < freq_specs[i + 1]) {
               scale = (freq[k] - freq_specs[i]) / (freq_specs[i + 1] - freq_specs[i]);

               for (int j = 0; j < 2; j++) {
                  gain[j][k] = 10 * log10(specs2[i][j] + scale * (specs2[i + 1][j] - specs2[i][j]));
               } //for
               i = 4;
            } //if
         } //for
      } //else
   } //for (frequencies)

   return 1;
}


double GaintoHeight(double gain, double freq)
{
   // from f=4*pi*A_eff/lambda^2
   // and h_eff=2*sqrt(A_eff*Z_rx/Z_air)
  //Shouldn't there be a factor 1/NFIRN^2 to shift frequencies?
   return 2.0*sqrt(gain / 4 / PI * C * C / (freq * freq) * Zr / Z0 * NICE);
} //GaintoHeight

double GetHeff(int AntType, double freq, double* n_boresight, double* n_epol, double* n_prop, double* n_pol)
{

  if (AntType==0){//Isotropic Antenne
    double Heff = GaintoHeight(gainv, freq * 1.E6);

    //    cout<<"AntType = "<<AntType<<", freq = "<<freq<<", Heff = "<<Heff<<endl;

    return Heff;
}
  else if (AntType==1){//Old SHelfMC LPDA

    double hitangle_e, hitangle_h, e_component, h_component;
    GetHitAngle(n_boresight, n_epol, n_prop, n_pol, hitangle_e, hitangle_h, e_component, h_component);

    double flare_i[4];
    GetFlare(freq, flare_i);

    double Heff = GaintoHeight(gainv, freq * 1.E6) *
           sqrt((pow(e_component * exp(-2 * ALOG2 * (hitangle_e / flare_i[0]) * (hitangle_e / flare_i[0])), 2)  +   pow(e_component * exp(-2 * ALOG2 * (hitangle_h / flare_i[1]) * (hitangle_h / flare_i[1])), 2)) / 2);

    //    cout<<"AntType = "<<AntType<<", freq = "<<freq<<", Heff = "<<Heff<<endl;

    return Heff;
  }
  else if (AntType==2){//100MHz Create LPDA with Anna's Antenna Model Framework

    int NC = Create100->N[0];

    freq = freq/1.30; //since this antenna model is in the firn, frequencies should be shifted

    double hitangle_e, hitangle_h, e_component, h_component;
    GetHitAngle(n_boresight, n_epol, n_prop, n_pol, hitangle_e, hitangle_h, e_component, h_component);//should use arrival direction, but it doesn't matter since we only use abs(e_component)

    double Re_Z = Create100->InterpolateToSingleFrequency(freq,  NC, Create100->frequencies, Create100->Re_Z);

    double gain = Create100->InterpolateToSingleFrequency(freq,  NC, Create100->frequencies, Create100->gains);

    /*
    cout<<"freq = "<<freq<<endl;
    cout<<"gain = "<<gain<<endl;
    cout<<"Re_Z = "<<Re_Z<<endl;
    cout<<"EffectiveHeight = "<<Create100->GetEffectiveHeight(gain,freq,Re_Z,C,119.99169*PI)<<endl;
    cout<<"abs(e_component) = "<<abs( e_component)<<endl;

*/
    //    double Heff = 2.0*C/(freq)*sqrt(gain / 4 / PI * C * C / (freq * freq) * Zr / Z0 * NICE)
    double Heff =  Create100->GetEffectiveHeight(gain,freq,Re_Z,C,119.99169*PI) * abs( e_component);

    //    cout<<"AntType = "<<AntType<<", freq = "<<freq<<", Heff = "<<Heff<<endl;

    return Heff;
  }
  else if (AntType==3){//ARA Bicone dipole

    int NC = ARA_Bicone->N[0];

    double Re_Z = ARA_Bicone->InterpolateToSingleFrequency(freq,  NC, ARA_Bicone->frequencies, ARA_Bicone->Re_Z);

    double gain = ARA_Bicone->InterpolateToSingleFrequency(freq,  NC, ARA_Bicone->frequencies, ARA_Bicone->gains);

    double Heff =  ARA_Bicone->GetEffectiveHeight(gain,freq,Re_Z,C,119.99169*PI) * abs(Dot(n_pol,n_epol));

    return Heff;

    }
    else if (AntType==4){//50MHz Create LPDA with Anna's Antenna Model Framework

      int NC = Create50->N[0];

      freq = freq/1.30; //since this antenna model is in the firn, frequencies should be shifted

      double hitangle_e, hitangle_h, e_component, h_component;
      GetHitAngle(n_boresight, n_epol, n_prop, n_pol, hitangle_e, hitangle_h, e_component, h_component);//should use arrival direction, but it doesn't matter since we only use abs(e_component)

      double Re_Z = Create50->InterpolateToSingleFrequency(freq,  NC, Create50->frequencies, Create50->Re_Z);

      double gain = Create50->InterpolateToSingleFrequency(freq,  NC, Create50->frequencies, Create50->gains);

      double Heff =  Create50->GetEffectiveHeight(gain,freq,Re_Z,C,119.99169*PI) * abs( e_component);

      return Heff;
    }
    else if (AntType==5){//100MHz Create LPDA In Firn with simulated boundary

      int NC = Create100Boundary->N[0];

      freq = freq/1.30; //since this antenna model is in the firn, frequencies should be shifted

      double hitangle_e, hitangle_h, e_component, h_component;
      GetHitAngle(n_boresight, n_epol, n_prop, n_pol, hitangle_e, hitangle_h, e_component, h_component);//should use arrival direction, but it doesn't matter since we only use abs(e_component)

      double Re_Z = Create100Boundary->InterpolateToSingleFrequency(freq,  NC, Create100Boundary->frequencies, Create100Boundary->Re_Z);

      double gain = Create100Boundary->InterpolateToSingleFrequency(freq,  NC, Create100Boundary->frequencies, Create100Boundary->gains);

      double Heff =  Create100Boundary->GetEffectiveHeight(gain,freq,Re_Z,C,119.99169*PI) * abs( e_component);

      return Heff;
    }
  else {
    cout<<"invalid Antenna Type!"<<endl;
    return -1;
  }
}

void GetPolarization(double* nnu, double* posnu2At, double* n_pol, double* n_Bfield)
{
   double Bfield[3];
   double Efield[3];
   Cross(nnu, posnu2At , Bfield);
   Cross(Bfield, posnu2At, Efield);

   nVector(Bfield, n_Bfield); //KD 07/19/11, our normalized Bfield
   nVector(Efield, n_pol);
   //check whether Efiled is normalized
   // cout<<Mag(nEfield)<<endl;
}


void GetRotation(double theta1, double theta2, double* original_nsignal, double* n_pol, double* n_pol_out)
{
//We use this technique for graded index, because in that case, there's no transmission loss
//only the polarization vector gets rotated smoothly

//Here, we use the 3D rotation matrix where the rotation axis is already defined as follows

   double rot_angle = theta2 - theta1;
   double rot_axis[3];
   double n_rot_axis[3];

//9/5/11; after a lot of troubleshooting, I figured out that rotation axis is flipped for reflected, ie downgoing signals.

   if (original_nsignal[2] > 0) {
      //upgoing direct signals have a different axis
      rot_axis[0] = -original_nsignal[1];
      rot_axis[1] = +original_nsignal[0];
   } else { // for reflected, there's a different rotation axis
      rot_axis[0] = +original_nsignal[1];
      rot_axis[1] = -original_nsignal[0];
   }
   rot_axis[2] = 0; //we make the axis lie on the x-y plane

   nVector(rot_axis, n_rot_axis); //need a normalized rotation axis

//now, evaluating each component of the 3x3 rotation matrix with the incident pol vector
   double pol_out[3]; // now defining our polarization in firn, presumably non-normalized
   pol_out[0] =   n_pol[0] * (n_rot_axis[0] * n_rot_axis[0] * (1 - cos(rot_angle)) + cos(rot_angle)) +
                  n_pol[1] * (n_rot_axis[0] * n_rot_axis[1] * (1 - cos(rot_angle)) - n_rot_axis[2] * sin(rot_angle)) +
                  n_pol[2] * (n_rot_axis[0] * n_rot_axis[2] * (1 - cos(rot_angle)) + n_rot_axis[1] * sin(rot_angle)) ;

   pol_out[1] =   n_pol[0] * (n_rot_axis[0] * n_rot_axis[1] * (1 - cos(rot_angle)) + n_rot_axis[2] * sin(rot_angle)) +
                  n_pol[1] * (n_rot_axis[1] * n_rot_axis[1] * (1 - cos(rot_angle)) + cos(rot_angle)) +
                  n_pol[2] * (n_rot_axis[1] * n_rot_axis[2] * (1 - cos(rot_angle)) - n_rot_axis[0] * sin(rot_angle)) ;

   pol_out[2] =   n_pol[0] * (n_rot_axis[0] * n_rot_axis[2] * (1 - cos(rot_angle)) - n_rot_axis[1] * sin(rot_angle)) +
                  n_pol[1] * (n_rot_axis[1] * n_rot_axis[2] * (1 - cos(rot_angle)) + n_rot_axis[0] * sin(rot_angle)) +
                  n_pol[2] * (n_rot_axis[2] * n_rot_axis[2] * (1 - cos(rot_angle)) + cos(rot_angle)) ;

   nVector(pol_out, n_pol_out);
}


//08/25/11 KD: this is adapted from shared.cc from ANITA MC
// I redefined Vectors to arrays to keep in line with main code
void GetFresnel(//const Vector &surface_normal,
   double theta1,//double *incident_angle,
   double theta2,//double *transmitted_angle,
   double* original_nsignal,//double *incident_signal[3], //const Vector &firn_rf,
   double* nsignal_atAT, //const Vector &air_rf
   double* n_pol, //Vector &pol,
   double* n_pol_out
   //double efield,
   //double &t_coeff_pokey,double &t_coeff_slappy,
   //double &fresnel,
   //double &mag
)
{


   // find angle of incidence and angle of transmission
//double incident_angle = surface_normal.Angle(firn_rf);
//double transmitted_angle = surface_normal.Angle(air_rf);
   double incident_angle = theta1;
   double transmitted_angle = theta2;

// cout<<">>angles "<<incident_angle<<" "<<transmitted_angle<<endl;

   double surface_normal[3] = {0, 0, 1}; //our regular positive z-axis
   double perp[3];
   double firn_parallel[3];
   double ice_parallel[3];
   double t_coeff_pokey, t_coeff_slappy;

//cout<<">>original_nsignal "<<original_nsignal[0]<<" "<<original_nsignal[1]<<" "<<orignal_nsignal[2]<<endl;
//cout<<">>nsignal_atAT "<<nsignal_atAT[0]<<" "<<nsignal_atAT[1]<<" "<<nsignal_atAT[2]<<endl;

   // this is perp the surface normal and transmitted ray, parallel to surface
   Cross(nsignal_atAT, surface_normal, perp);
   //Vector perp = air_rf.Cross(surface_normal).Unit();
//cout<<">>perp "<<perp[0]<<" "<<perp[1]<<" "<<perp[2]<<endl;

   // this is in the bending plane
   Cross(perp, nsignal_atAT, firn_parallel);
   //Vector air_parallel = perp.Cross(air_rf).Unit();
//cout<<">>firn_parallel "<<firn_parallel[0]<<" "<<firn_parallel[1]<<" "<<firn_parallel[2]<<endl;

   // this is in the bending plane
   Cross(perp, original_nsignal, ice_parallel); //not air->firn, and firn->ice
   //Vector firn_parallel = perp.Cross(firn_rf).Unit();
//cout<<">>ice_parallel "<<ice_parallel[0]<<" "<<ice_parallel[1]<<" "<<ice_parallel[2]<<endl;

   // component of polarization (in the firn) perp to surface normal
   double pol_perp_ice = Dot(n_pol, perp); // this is the slappy component in the ice
   double pol_parallel_ice = Dot(n_pol, ice_parallel); // this is the pokey component in the ice
//cout<<">>pol_perp_ice "<<pol_perp_ice<<" pol_parallel_ice "<<pol_parallel_ice<<endl;

   double pol_perp_firn, pol_parallel_firn;
   double r_coeff_pokey =  tan(incident_angle - transmitted_angle) / tan(incident_angle + transmitted_angle);
   t_coeff_pokey = sqrt((1. - r_coeff_pokey * r_coeff_pokey));
   pol_parallel_firn = t_coeff_pokey * pol_parallel_ice; // find pokey component in the firn/air

   double r_coeff_slappy = sin(incident_angle - transmitted_angle) / sin(incident_angle + transmitted_angle);
   t_coeff_slappy = sqrt((1. - r_coeff_slappy * r_coeff_slappy));
   pol_perp_firn = t_coeff_slappy * pol_perp_ice; // find slappy component in the ice/firn
//cout<<">>r_coeffs: "<<r_coeff_pokey<<" "<<r_coeff_slappy<<endl;
//cout<<">>pol_perp_firn "<<pol_perp_firn<<" pol_parallel_firn "<<pol_parallel_firn<<endl;


//THIS CAN BE USED FOR MODIFYING THE AMPLITUDES
   //mag=sqrt( tan(incident_angle) / tan(transmitted_angle) );
   //fresnel = sqrt( pow(efield * pol_perp_firn,2) + pow(efield * pol_parallel_firn,2)) / efield;


   double pol_out[3]; //double n_pol_out[3];
//pol_out = (pol_perp_firn * perp + pol_parallel_firn * firn_parallel);//.Unit();

   pol_out[0] = (pol_perp_firn * perp[0] + pol_parallel_firn * firn_parallel[0]);//.Unit();
   pol_out[1] = (pol_perp_firn * perp[1] + pol_parallel_firn * firn_parallel[1]);
   pol_out[2] = (pol_perp_firn * perp[2] + pol_parallel_firn * firn_parallel[2]);


//cout<<">>n_pol "<<n_pol[0]<<" "<<n_pol[1]<<" "<<n_pol[2]<<endl;
//cout<<">>pol_out "<<pol_out[0]<<" "<<pol_out[1]<<" "<<pol_out[2]<<endl;
   nVector(pol_out, n_pol_out);
//cout<<">>n_pol_out "<<n_pol_out[0]<<" "<<n_pol_out[1]<<" "<<n_pol_out[2]<<endl;

} //GetFresnel


void GetHitAngle(double* n_boresight, double* n_epol, double* n_prop, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component,
                     double& h_component)
{

   double n_hplane[3];
   Cross(n_boresight,n_epol,n_hplane);

   //these are the components of the propagation vector in the e and h plane.  They are overwritten later in the function
   e_component = -Dot(n_prop, n_epol);
   h_component = -Dot(n_prop, n_hplane);
   double n_component = -Dot(n_prop, n_boresight);

   hitangle_e = atan(h_component / n_component);
   if (n_component < 0 && h_component < 0)
      hitangle_e -= PI;
   if (n_component < 0 && h_component > 0)
      hitangle_e += PI;

   hitangle_h = atan(e_component / n_component);
   if (n_component < 0 && e_component < 0)
      hitangle_h -= PI;
   if (n_component < 0 && e_component > 0)
      hitangle_h += PI;

   //NOTE: this is a re-evaluated e_component and h_component that goes into term_LPA calculation
   e_component = Dot(n_pol, n_epol);
   h_component = Dot(n_pol, n_hplane);

}


void GetHitAngle_ST0(double* nposnu2AT, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component,
                     double& h_component)
{


   double n_normal[3] = {0, 0, -1};
   double n_hplane[3] = {1, 0, 0};
   double n_epol[3] = {0, 1, 0};
   double n_component;
   e_component = -Dot(nposnu2AT, n_epol);
   h_component = -Dot(nposnu2AT, n_hplane);
   n_component = -Dot(nposnu2AT, n_normal);
   hitangle_e = atan(h_component / n_component);
   if (n_component < 0 && h_component < 0)
      hitangle_e -= PI;
   if (n_component < 0 && h_component > 0)
      hitangle_e += PI;
   hitangle_h = atan(e_component / n_component);
   if (n_component < 0 && e_component < 0)
      hitangle_h -= PI;
   if (n_component < 0 && e_component > 0)
      hitangle_h += PI;
   e_component = Dot(n_pol, n_epol);
   h_component = Dot(n_pol, n_hplane);

}


void GetHitAngle_ST2(double* nposnu2AT, double* n_pol, double* hitangle_e, double* hitangle_h, double* e_component,
                     double* h_component)
//KD: Z is downward angle specified in input file
{

   double  n_normal[5][3], n_hplane[5][3], n_epol[5][3];
   double tempn[5][3] = {{0, 0, -1}, {cos(Z), 0, -sin(Z)}, {0, -cos(Z), -sin(Z)}, { -cos(Z), 0, -sin(Z)}, {0, cos(Z), -sin(Z)}};
   double temph[5][3] = {{1, 0, 0}, {0, -1, 0}, { -1, 0, 0}, {0, 1, 0}, {1, 0, 0}};
   double tempe[5][3] = {{0, 1, 0}, {sin(Z), 0, cos(Z)}, {0, -sin(Z), cos(Z)}, { -sin(Z), 0, cos(Z)}, {0, sin(Z), cos(Z)}};
   for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++) {
         n_normal[i][j] = tempn[i][j];
         n_hplane[i][j] = temph[i][j];
         n_epol[i][j] = tempe[i][j];
      }
   }



   double n_component[5];

   for (int i = 0; i < 5; i++) {
      e_component[i] = -Dot(nposnu2AT, n_epol[i]); //use the antenna as the tail of a vector
      h_component[i] = -Dot(nposnu2AT, n_hplane[i]);
      n_component[i] = -Dot(nposnu2AT, n_normal[i]);
      hitangle_e[i] = atan(h_component[i] / n_component[i]);
      if (n_component[i] < 0 && h_component[i] < 0)
         hitangle_e[i] -= PI;
      if (n_component[i] < 0 && h_component[i] > 0)
         hitangle_e[i] += PI;

      hitangle_h[i] = atan(e_component[i] / n_component[i]);
      if (n_component[i] < 0 && e_component[i] < 0)
         hitangle_h[i] -= PI;
      if (n_component[i] < 0 && e_component[i] > 0)
         hitangle_h[i] += PI;
      e_component[i] = Dot(n_pol, n_epol[i]);
      h_component[i] = Dot(n_pol, n_hplane[i]);
   }

}

void GetHitAngle_ST1(double* nposnu2AT, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component,
                     double& h_component, double& hitangle_e2, double& hitangle_h2, double& e2_component, double& h2_component)
{
   double n_epol[3] = {0, 1, 0}; //y direction,unit vector in e plane
   double n_hplane[3] = {1, 0, 0}; //x direction, unit vector in h plane
   double n_normal[3] = {0, 0, -1}; //-z direction, normal direction which is perpendicular to the face of an antenna
   double n_epol2[3] = {0.707, 0.707, 0};
   double n_hplane2[3] = {0.707, -0.707, 0};
   double n_normal2[3] = {0, 0, -1};
   //double mag_posnu2AT=Mag(posnu2AT);
   double n_component = 0;
   double n2_component = 0;
   e_component = -Dot(nposnu2AT, n_epol);
   e2_component = -Dot(nposnu2AT, n_epol2);
   h_component = -Dot(nposnu2AT, n_hplane);
   h2_component = -Dot(nposnu2AT, n_hplane2);
   n_component = -Dot(nposnu2AT, n_normal);
   n2_component = -Dot(nposnu2AT, n_normal2);


   hitangle_e = atan(h_component / n_component);
   if (n_component < 0 && h_component < 0)
      hitangle_e -= PI;
   if (n_component < 0 && h_component > 0)
      hitangle_e += PI;

   hitangle_e2 = atan(h2_component / n2_component);
   if (n2_component < 0 && h2_component < 0)
      hitangle_e2 -= PI;
   if (n2_component < 0 && h2_component > 0)
      hitangle_e2 += PI;

   hitangle_h = atan(e_component / n_component);
   if (n_component < 0 && e_component < 0)
      hitangle_h -= PI;
   if (n_component < 0 && e_component > 0)
      hitangle_h += PI;

   hitangle_h2 = atan(e2_component / n2_component);
   if (n2_component < 0 && e2_component < 0)
      hitangle_h2 -= PI;
   if (n2_component < 0 && e2_component > 0)
      hitangle_h2 += PI;

   e2_component = Dot(n_pol, n_epol2);
   h2_component = Dot(n_pol, n_hplane2);
   e_component = Dot(n_pol, n_epol);
   h_component = Dot(n_pol, n_hplane);


}

void GetHitAngle_LPA(int NofAT, int N_Ant_perST, double* nposnu2AT, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component, double& h_component) //KD: added N_Ant_perST
{
   //hitangle_e is the angle between the signal and the E-plane of LPA( which is the plane of the dipole array)
   //hitangle_h is the angle between the signal and the H-plane of LPA(which is perpendicular to the dipole array plane)
   //pol_component is the polarization component in the dipole's polarization direction
   double n_normal[3] = {0, 0, -1}; //vertically downward
   double n_epol[3] = {cos((0.5 + NofAT * (2. / N_Ant_perST))*PI), sin((0.5 + NofAT * (2. / N_Ant_perST))*PI), 0.};
   //FW n_epol[3]={cos((2.+NofAT)*0.25*PI),sin((2.+NofAT)*0.25*PI),0.};
   //the unit vector in the e-plane(also called the array plane)
   double n_hplane[3] = {cos((2. / N_Ant_perST)*PI * NofAT), sin((2. / N_Ant_perST)*PI * NofAT), 0.};
   //FW       n_hplane[3]={cos(0.25*PI*NofAT),sin(0.25*PI*NofAT),0.};
   //perpendicular to the array plane,point outward from the center of the station
   e_component = -Dot(nposnu2AT, n_epol);
   h_component = -Dot(nposnu2AT, n_hplane);
   double n_component = -Dot(nposnu2AT, n_normal);

//cout<<"KDAa: "<<"e,h,n components:"<<"["<<e_component<<" , "<<h_component<<" , "<<n_component<<"]"<<endl;

   hitangle_e = atan(h_component / n_component);
   if (n_component < 0 && h_component < 0)
      hitangle_e -= PI;
   if (n_component < 0 && h_component > 0)
      hitangle_e += PI;

   hitangle_h = atan(e_component / n_component);
   if (n_component < 0 && e_component < 0)
      hitangle_h -= PI;
   if (n_component < 0 && e_component > 0)
      hitangle_h += PI;

   //NOTE: this is a re-evaluated e_component and h_component that goes into term_LPA calculation
   e_component = Dot(n_pol, n_epol);
   h_component = Dot(n_pol, n_hplane);

//cout<<"KDAb: "<<"e,h,n components:"<<"["<<e_component<<" , "<<h_component<<" , "<<n_component<<"]"<<endl;

}

void GetMirrorHitAngle_LPA(int NofAT, int N_Ant_perST, double* nposnu2MirrorAT, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component, double& h_component)
{
   //hitangle_e is the angle between the signal and the E-plane of LPA( which is the plane of the dipole array)
   //hitangle_h is the angle between the signal and the H-plane of LPA(which is perpendicular to the dipole array plane)
   //pol_component is the polarization component in the dipole's polarization direction
   double n_normal[3] = {0, 0, 1}; //vertically upward
//FW------------------------------------------------------
//  double n_epol[3]={cos((2.+NofAT)*0.25*PI),sin((2.+NofAT)*0.25*PI),0.};
   //the unit vector in the e-plane(also called the array plane)
//  double n_hplane[3]={cos(0.25*PI*NofAT),sin(0.25*PI*NofAT),0.};
   //perpendicular to the array plane,point outward from the center of the station
//--------------------------------------------------------
   double n_epol[3] = {cos((0.5 + NofAT * (2. / N_Ant_perST))*PI), sin((0.5 + NofAT * (2. / N_Ant_perST))*PI), 0.};
   //the unit vector in the e-plane(also called the array plane)
   double n_hplane[3] = {cos((2. / N_Ant_perST)*PI * NofAT), sin((2. / N_Ant_perST)*PI * NofAT), 0.};
   //perpendicular to the array plane,point outward from the center of the station

   e_component = -Dot(nposnu2MirrorAT, n_epol);
   h_component = -Dot(nposnu2MirrorAT, n_hplane);
   double n_component = -Dot(nposnu2MirrorAT, n_normal);
   hitangle_e = atan(h_component / n_component);
   if (n_component < 0 && h_component < 0)
      hitangle_e -= PI;
   if (n_component < 0 && h_component > 0)
      hitangle_e += PI;
   hitangle_h = atan(e_component / n_component);
   if (n_component < 0 && e_component < 0)
      hitangle_h -= PI;
   if (n_component < 0 && e_component > 0)
      hitangle_h += PI;

   e_component = Dot(n_pol, n_epol);
   h_component = Dot(n_pol, n_hplane);
}

void GetMirrorHitAngle_ST0(double* n_posnu2MirrorAT, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component, double& h_component)
{
   double n_normal[3] = {0, 0, 1};
   double n_hplane[3] = {1, 0, 0};
   double n_epol[3] = {0, 1, 0};
   double n_component;
   e_component = -Dot(n_posnu2MirrorAT, n_epol);
   h_component = -Dot(n_posnu2MirrorAT, n_hplane);
   n_component = -Dot(n_posnu2MirrorAT, n_normal);
   hitangle_e = atan(h_component / n_component);
   if (n_component < 0 && h_component < 0)
      hitangle_e -= PI;
   if (n_component < 0 && h_component > 0)
      hitangle_e += PI;
   hitangle_h = atan(e_component / n_component);
   if (n_component < 0 && e_component < 0)
      hitangle_h -= PI;
   if (n_component < 0 && e_component > 0)
      hitangle_h += PI;
   e_component = Dot(n_pol, n_epol);
   h_component = Dot(n_pol, n_hplane);
}

void GetMirrorHitAngle_ST2(double* n_posnu2MirrorAT, double* n_pol, double* hitangle_e, double* hitangle_h, double* e_component, double* h_component)
{

   double n_normal[5][3], n_hplane[5][3], n_epol[5][3];
   double tempn[5][3] = {{0, 0, 1}, {cos(Z), 0, sin(Z)}, {0, -cos(Z), sin(Z)}, { -cos(Z), 0, sin(Z)}, {0, cos(Z), sin(Z)}}; //all the z components of n turn to negative of the original
   double temph[5][3] = {{1, 0, 0}, {0, -1, 0}, { -1, 0, 0}, {0, 1, 0}, {1, 0, 0}}; //horizontal component doesn't change
   double tempe[5][3] = {{0, 1, 0}, {sin(Z), 0, -cos(Z)}, {0, -sin(Z), -cos(Z)}, { -sin(Z), 0, -cos(Z)}, {0, sin(Z), -cos(Z)}}; //change z components to be -z
   for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 3; j++) {
         n_normal[i][j] = tempn[i][j];
         n_hplane[i][j] = temph[i][j];
         n_epol[i][j] = tempe[i][j];
      }
   }


   double n_component[5];

   for (int i = 0; i < 5; i++) {
      e_component[i] = -Dot(n_posnu2MirrorAT, n_epol[i]); //use the antenna as the tail of a vector
      h_component[i] = -Dot(n_posnu2MirrorAT, n_hplane[i]);
      n_component[i] = -Dot(n_posnu2MirrorAT, n_normal[i]);
      hitangle_e[i] = atan(h_component[i] / n_component[i]);
      if (n_component[i] < 0 && h_component[i] < 0)
         hitangle_e[i] -= PI;
      if (n_component[i] < 0 && h_component[i] > 0)
         hitangle_e[i] += PI;

      hitangle_h[i] = atan(e_component[i] / n_component[i]);
      if (n_component[i] < 0 && e_component[i] < 0)
         hitangle_h[i] -= PI;
      if (n_component[i] < 0 && e_component[i] > 0)
         hitangle_h[i] += PI;

      //  cout<<hitangle_e[i]*180/PI<<"  mirror       eeeeeeeeeeeee"<<endl;
      // cout<<hitangle_h[i]*180/PI<<"  mirror       hhhhhhhhhhhhh"<<endl;

      e_component[i] = Dot(n_pol, n_epol[i]);
      h_component[i] = Dot(n_pol, n_hplane[i]);
   }



}

void GetMirrorHitAngle_ST1(double* nposnu2MirrorAT, double* n_pol, double& hitangle_e, double& hitangle_h, double& e_component, double& h_component, double& hitangle_e2, double& hitangle_h2, double& e2_component, double& h2_component)
{
   double n_epol[3] = {0, 1, 0}; //y direction,unit vector in e plane
   double n_hplane[3] = {1, 0, 0}; //x direction, unit vector in h plane
   double n_normal[3] = {0, 0, 1}; //-z direction, normal direction which is perpendicular to the face of an antenna
   double n_epol2[3] = {0.707, 0.707, 0};
   double n_hplane2[3] = {0.707, -0.707, 0};
   double n_normal2[3] = {0, 0, 1};
   double n_component = 0;
   double n2_component = 0;
   e_component = -Dot(nposnu2MirrorAT, n_epol);
   e2_component = -Dot(nposnu2MirrorAT, n_epol2);
   h_component = -Dot(nposnu2MirrorAT, n_hplane);
   h2_component = -Dot(nposnu2MirrorAT, n_hplane2);
   n_component = -Dot(nposnu2MirrorAT, n_normal);
   n2_component = -Dot(nposnu2MirrorAT, n_normal2);

   hitangle_e = atan(h_component / n_component);
   if (n_component < 0 && h_component < 0)
      hitangle_e -= PI;
   if (n_component < 0 && h_component > 0)
      hitangle_e += PI;

   hitangle_e2 = atan(h2_component / n2_component);
   if (n2_component < 0 && h2_component < 0)
      hitangle_e2 -= PI;
   if (n2_component < 0 && h2_component > 0)
      hitangle_e2 += PI;

   hitangle_h = atan(e_component / n_component);
   if (n_component < 0 && e_component < 0)
      hitangle_h -= PI;
   if (n_component < 0 && e_component > 0)
      hitangle_h += PI;

   hitangle_h2 = atan(e2_component / n2_component);
   if (n2_component < 0 && e2_component < 0)
      hitangle_h2 -= PI;
   if (n2_component < 0 && e2_component > 0)
      hitangle_h2 += PI;



   // hitangle_e=fabs(PI/2-acos(e_component/Mag(posnu2AT)/Mag(n_epol)));
   // hitangle_e2=fabs(PI/2-acos(e2_component/Mag(posnu2AT)/Mag(n_epol2)));
   // hitangle_h=fabs(PI/2-acos(h_component/Mag(posnu2AT)/Mag(n_hplane)));
   // hitangle_h2=fabs(PI/2-acos(h2_component/Mag(posnu2AT)/Mag(n_hplane2)));
   e2_component = Dot(n_pol, n_epol2);
   h2_component = Dot(n_pol, n_hplane2);
   e_component = Dot(n_pol, n_epol);
   h_component = Dot(n_pol, n_hplane);


}

void ConvertEHtoLR(double e_component, double h_component, double& lcp_component, double& rcp_component)
{

   lcp_component = sqrt((e_component * e_component + h_component * h_component) / 2);
   rcp_component = sqrt((e_component * e_component + h_component * h_component) / 2);

} //ConvertEHtoLR

int Getifreq(double freq)
{

   return (int)((freq - FREQ_LOW) / (FREQ_HIGH - FREQ_LOW) * (double)NFREQ);
} //Getifreq



double GetN(double height)
//KD 07/21/11 I need to rewrite this so that when FIRN is 0 in input file, I don't have to systematically turn off DEPTH_DEPENDENT in declaration.hh file
{
    double n = 0;

    if (!DEPTH_DEPENDENT_N) {
       if (FIRN && height > (ICETHICK - FIRNDEPTH))
          n = NFIRN;
       else
          n = NICE;
    }
    else {
       //these are Peter's fit parameters//old usage by FW
       //double a1=0.463251;
       //double b1=0.0140157;
       //n=NFIRN+a1*(1.0-exp(-b1*(ICETHICK-posnu[2]))); //KD: this was the previous one

       //this is from BARELLA paper, use this now
       //n = 1.0 + 0.86*(1.0-0.638*exp(-(ICETHICK-posnu[2])/34.7));

       //n = 1.0 + 0.86 * (1.0 - 0.638 * exp(-FIRNfactor * (ICETHICK - posnu[2]) / 34.7));
       n = NICE - (NICE - NFIRN) * exp(-(ICETHICK - height) / C_INDEX);
     }
   return n;

}

void CalcShadowEdge(double zstep)
{
    double SinTheta0 = NFIRN/NICE;
    double SinSquareTheta0 = SinTheta0*SinTheta0;
    double Nz = NICE;
    double dr = 0;
    double rtemp = 0;
    double z;
    z = zstep;

    SHADOWRANGE.push_back(0.0);
    SHADOWHEIGHTS.push_back(0.0);

    while (z <= ICETHICK)
    {
        Nz = GetN(z);
        if (Nz>NFIRN){
          dr = zstep*SinTheta0/sqrt(Nz*Nz/NICE/NICE - SinSquareTheta0);
        }
        else{
          dr=58*zstep; //treat a horizontal ray like 1deg
        }
        if (dr > 58*zstep){
          //The idea is to keep the ray steeper than 1degree, keeping 1/tan from blowing up
          rtemp += 58*zstep;
        }
        else{
          rtemp += dr;
        }
        SHADOWRANGE.push_back(rtemp);
        SHADOWHEIGHTS.push_back(z);
        z+=zstep;
    }

    int nPoints = SHADOWRANGE.size();

    for (int i=0;i<nPoints;i++)
    {
        SHADOWRANGE[i] = rtemp - SHADOWRANGE[i];
    }
    if (SHADOWHEIGHTS.back()<ICETHICK){
      SHADOWRANGE.push_back(0.0);
      SHADOWHEIGHTS.push_back(ICETHICK);
    }
}

double InterpolateLinear(double xtest, double x0, double y0, double x1, double y1)
//must have x0<=xtest<=x1, x0!=x1
{
    return y0 + (y1-y0)/(x1-x0)*(xtest-x0);
}

double GetRange(double height)
{
    //cout << "\nGetRange for height = " << height << endl;
    int npoints = SHADOWHEIGHTS.size();
    int i = int(npoints*height/ICETHICK);
    //cout << "initial guess i = " << i << endl;
    double heighti = SHADOWHEIGHTS[i];
    //cout << "initial guess height = " << heighti << endl;
    double rangei = SHADOWRANGE[i];
    //cout << "initial guess range = " << rangei << endl;
    if (height == heighti){
      //cout << "correct height guess, return range = " << rangei << endl << endl;
      return rangei;
    }
    else {
      int lasti = i;
      double heightlasti = SHADOWHEIGHTS[lasti];
      bool greater = height > heighti;
      //cout << "greater = " << greater << endl;
      bool changed = false;
      while (changed == false){
          lasti = i;
          heightlasti = heighti;
          if (greater)
            i++;
          else
            i--;
          if (i>=npoints){
            return 0.0;
          }
          else{
            heighti=SHADOWHEIGHTS[i];
            //cout << "lasti = " << lasti << ", i = " << i << ", heightlasti = " << heightlasti << ", heighti = " << heighti << endl;
            if ((heighti - height)*(heightlasti - height) <= 0){
                changed = true;
              }
          }
      }
      double rangelasti = SHADOWRANGE[lasti];
      rangei = SHADOWRANGE[i];
      //cout << "rangelasti = " << rangelasti << "rangei = " << rangei << endl << endl;
      if (greater){
          return InterpolateLinear(height,heighti,rangei,heightlasti,rangelasti);
      }
      else{
          return InterpolateLinear(height,heightlasti,rangelasti,heighti,rangei);
      }
    }
}

/*//Now calculated for a custom index profile
double GetRange(double height)
//KD: 9/23/2010 introducing a variable horizontal range due to shadowing
//KD need to have this activated only if FIRN present
{
    double depth = ICETHICK - height;
    if (depth<=FIRNDEPTH){
        return 13.09 * pow(depth, 0.56) - 0.74; //fit to range integration
    }
    else{
        return 13.09 * pow(FIRNDEPTH, 0.56) - 0.74 + (depth-FIRNDEPTH); //below firn, extend at 45 degrees
    }

}
*/

void Zero(double* array, int i)
{
   for (int k = 0; k < i; k++) {
      array[k] = 0;
   }

}
void Zero(int* array, int i)
{
   for (int k = 0; k < i; k++) {
      array[k] = 0;
   }
}
double fx(double x, double h1, double h2, double nconst, double deltax)
{
   return h1 / sqrt(1 / x - 1) + h2 / sqrt(nconst / x - 1) - deltax;
}

void Refract(double deltax, double yIce, double yFirn, double nVertex, double nSurface, double &dIce, double &dFirn, double &theta1, double &theta2)
{
  double x1 = 1.e-100;
  double x3 = 0.;
  double nconst = nSurface * nSurface / nVertex / nVertex; //here nSurface=NFIRN, nVertex=NICE
  double x2 = nconst; // x2 is not angle. x2=sin(theta)^2

  if (deltax == 0) {
     theta1 = 0.;
     theta2 = 0.;
     dIce = yIce;
     dFirn = yFirn;
  }

  else {
     do {
        x3 = (x1 + x2) / 2.;
        if (fx(x1, yIce, yFirn, nconst, deltax)*fx(x3, yIce, yFirn, nconst, deltax) < 0)
           x2 = x3;
        else x1 = x3;
     }while (fabs(fx(x3, yIce, yFirn, nconst, deltax)) > 0.00001);

     theta1 = asin(sqrt(x3));
     theta2 = asin(nVertex * sin(theta1) / nSurface);
     dIce = yIce / cos(theta1);
     dFirn = yFirn / cos(theta2);
  }

}

double GetVmMHz1m(double pnu, double freq, double X0DEPTH, double ECDEPTH, double NDEPTH, double AEXDEPTH)
{

   // parametrization from Jaime Alvarez Munhiz
   //  here using astro-ph/0003315
   double nu0 = 1150.;
   double nu0_modified = nu0
                         * (X0ICE / ECICE) / (X0DEPTH / ECDEPTH)
                         * (1 / sqrt(NDEPTH * NDEPTH - 1.)) / (1 / sqrt(NICE * NICE - 1.));

   // frequency in MHz

   double factor =
      //1/sin(changle) // should be cerenkov angle for ice
      1 / sqrt(1 - 1 / (NICE * NICE)) // sin(changle) for ice
      * 1 / nu0 //
      *X0DEPTH / X0ICE // track length *** use dE/dX rho instead
      *ECICE / ECDEPTH
      * AEXDEPTH / AEX_ICE; // to account for critical energy
   // to account for cerenkov threshold // maybe should be "a" instead

   double vmmhz1m = factor * (2.53E-7) * (pnu / 1.E12) * freq * (1. / (1. + pow(freq / nu0_modified, 1.44))) * JAIME_FACTOR;


   vmmhz1m /= sqrt(2.); //from ANITA's decision
   return vmmhz1m;
} //GetVmMHz1m

void GetVmMHz(double vmmhz_max, double vmmhz1m_max, double pnu, double* vmmhz)
{

   // parametrization from Jaime Alvarez Munhiz
   //  here using astro-ph/0003315

   for (int i = 0; i < NFREQ; i++) {

      vmmhz[i] = vmmhz_max / vmmhz1m_max * GetVmMHz1m(pnu, freq[i], X0ICE, ECICE, n1, AEX_ICE);

   }
} //GetVmMHz

void GetVmMHz_freq(double* vmmhz_max_freq, double vmmhz1m_max, double pnu, double* vmmhz)
{
   for (int i = 0; i < NFREQ; i++) {

      vmmhz[i] = vmmhz_max_freq[i] / vmmhz1m_max * GetVmMHz1m(pnu, freq[i], X0ICE, ECICE, n1, AEX_ICE);

   }

}

void GetSpread(double pnu,
               double emfrac,
               double hadfrac,
               double freq,
               double n_depth,
               double X0DEPTH,

               double& deltheta_em_max,
               double& deltheta_had_max,

               double& eshower_em, //KD
               double& eshower_had, //KD
               double& showerlength_em //KD
              )
{

   //  scale by how far off Cherenkov angle this viewing antenna is
   //  c.f. A-MZ  astro-ph/9706064 and astro-ph/0003315
   //  and for non-LPM (non-EM) showers from Phys.Lett.B434,396 (1998)
   //  The lengths are different hence the angular thickness of
   //  the shower is different.  Get the angular thickness for
   //  both the EM and hadroic parts.

   double elpm = GetLPM(X0DEPTH);

   //  freq=freq/1.E6;  // frequency in MHz
   double showerlength = 3.1; //shower length in meters-gets a modification
   //for em showers due to lpm effect.
   // this shower length is chosen somewhat arbitrarily, but is
   // approximately the length of a shower in ice.
   // Then, the coefficient out front of the equations for
   // deltheta_em_max and deltheta_had_max are set so that
   // for ice, we get the equations in astro-ph/9706064
   // with the coefficient in front being 2.7 degrees.
   // I wanted to make the dependence on the shower length
   // and index of refraction explicit, so I pulled those variables
   // out of the equations for deltheta_em_max and deltheta_had_max.

   double eshower;  // shower energy

   double nu0 = 500.E6 / 1.E6 * X0ICE / X0DEPTH; // for rego (astro-ph/9706064)
   // decoherence frequency scales with radiation length
   // when X0DEPTH=X0ICE, nu0=500 MHz as in astro-ph/9706064

   eshower = emfrac * pnu; // first, consider the electromagnetic shower.
   eshower_em = eshower;  //KD



   // lengthen the shower to account for the lpm effect.
   if (eshower < elpm || !LPM)
      showerlength /= pow((eshower / elpm), -0.03);
   else
      showerlength /= pow(elpm / (0.14 * (eshower) + elpm), 0.3);

   // these equations are in astro-ph/9706064, but we have pulled
   // out the dependence on index of refraction and shower length.
   // note that 12.32/sqrt(pow(n_depth,2)-1)*DEG2RAD/showerlength=2.7 degrees.

   showerlength_em = showerlength;

   deltheta_em_max = 12.32 / sqrt(pow(n_depth, 2) - 1) * (nu0 / freq) * DEG2RAD / showerlength;

//KD:02/16/11 these two lines were in the loop below, but it wasn't calculating the right eshower_had if y<0.00001, so I placed it outside now.
   eshower = hadfrac * pnu;
   eshower_had = eshower;

   if (hadfrac > 0.00001) { // if there is a hadronic component

      //eshower=hadfrac*pnu;  // just the energy of the hadronic component of the shower
//eshower_had = eshower;  //KD

      // these equations are in astro-ph/9706064, but we have pulled
      // out the dependence on index of refraction and shower length.
      double epsilon = log10(eshower / 1.E12);
      if (eshower >= 1E12 && eshower < 100.E12)
         deltheta_had_max = 1.473 / sqrt(pow(n_depth, 2) - 1) * nu0 / freq * DEG2RAD * (2.07 - 0.33 * epsilon + (7.5e-2) * epsilon * epsilon);
      else if (eshower < 100.E15)
         deltheta_had_max = 1.473 / sqrt(pow(n_depth, 2) - 1) * nu0 / freq * DEG2RAD * (1.744 - (1.21e-2) * epsilon);
      else if (eshower < 10.E18)
         deltheta_had_max = 1.473 / sqrt(pow(n_depth, 2) - 1) * nu0 / freq * DEG2RAD * (4.23 - 0.785 * epsilon + (5.5e-2) * epsilon * epsilon);
      else {
         //  beyond param, just use value at 10 EeV since slow variation
         //  and parameterization might diverge
         //  so scale from 10 EeV at 7.5% per decade (30/4=7.5)
         deltheta_had_max = 1.473 / sqrt(pow(n_depth, 2) - 1) * nu0 / freq * DEG2RAD * (4.23 - 0.785 * 7. + 5.5e-2 * 49.); // the last part in parenthesis if the previous equation evaluated at epsilon=7.
         deltheta_had_max = deltheta_had_max * (1. + (epsilon - 7.) * 0.075);
      } //else : beyond paramatrization
   } //if (hadronic component)

   else
      deltheta_had_max = 1.E-10;

}

void TaperVmMHz(double viewangle,
                double deltheta_em,
                double deltheta_had,
                double emfrac,
                double hadfrac,

                double& vmmhz,
                double& vmmhz_em, //KD for energy tracking
                double& vmmhz_had //KD for energy tracking
               )
{
   //--EM
   double delangle2 = (viewangle - changle) * (viewangle - changle);

//double vmmhz_em;
//double vmmhz_had;

   double rtemp = 0.5 * delangle2 / (deltheta_em * deltheta_em); //KD 01/21/11 should this be 0.5 or ln2?

   if (emfrac != 0) {
      if (rtemp <= 20) {
         vmmhz_em = vmmhz * exp(-rtemp);
         //if (fabs(viewangle-CHANGLE)/deltheta_em>2.7)
         //vmmhz_em=vmmhz/1000;
         //if (fabs(viewangle-CHANGLE)>0.44)
         //vmmhz_em=0;
      } else
         vmmhz_em = 0.;
   } else
      vmmhz_em = 0;


   //--HAD
   rtemp = ALOG2 * delangle2 / (deltheta_had * deltheta_had);

   if (hadfrac != 0) {
      if (rtemp < 20) {
         vmmhz_had = vmmhz * exp(-rtemp);
         //if (fabs(viewangle-CHANGLE)/deltheta_had>2.7)
         //vmmhz_had=vmmhz/1000;
         //if (fabs(viewangle-CHANGLE)>0.44)
         //vmmhz_had=0;
      } else
         vmmhz_had = 0.;
   } else
      vmmhz_had = 0.;

   vmmhz = sin(viewangle) * (emfrac * vmmhz_em + hadfrac * vmmhz_had);

}


void GetE2()
{
   double E[] = { -13.7, -14.2, -14.7, -15.2, -15.7, -16.2, -16.7, -17.2, -17.7, -18.2, -18.7, -19.2}; //from 16 to 21.5

   for (int i = 0; i < 12; i++) {
      energy[i] = 16. + ((double)i) / 2.;
   }
   for (int i = 0; i < 12; i++) {
      EdNdEdAdt[i] = pow(10., E[i]);
   }
   maxflux = Max(EdNdEdAdt, 12);
   for (int i = 0; i < 12; i++) {
      EdNdEdAdt[i] = EdNdEdAdt[i] / maxflux;
   }

}

void GetGZK()
{
   double Emuons[12];// E^2 dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
   //KD I think it's EdN/dE/dA/dt
   double Eelectrons[12];// E^2 dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
   //KD I think it's EdN/dE/dA/dt

   for (int i = 0; i < 12; i++) {
      energy[i] = 16. + ((double)i) / 2.;
      Emuons[i] = -30.;
      Eelectrons[i] = -30.;
   }


   /*
   // astro-ph.HE/1011.6340v2
   //Using fig 3(f), setting Eelectrons as representative of all flavors, and punting where needed
   //Note this is for point source, so we have to limit isotropic flux over 2pi stereadian
   Eelectrons[0]=log10(pow(10,-17.0)/2./PI); // 16.
   Eelectrons[1]=log10(pow(10,-17.186)/2./PI); // 16.5
   Eelectrons[2]=log10(pow(10,-17.371)/2./PI); // 17.
   Eelectrons[3]=log10(pow(10,-17.5)/2./PI); // 17.5
   Eelectrons[4]=log10(pow(10,-17.571)/2./PI); // 18.
   Eelectrons[5]=log10(pow(10,-18.071)/2./PI); // 18.5
   Eelectrons[6]=log10(pow(10,-19.571)/2./PI); // 19
   Eelectrons[7]=log10(pow(10,-24.9)/2./PI); // 19.5
   Eelectrons[8]=log10(pow(10,-30.)/2./PI); // 20.
   Eelectrons[9]=log10(pow(10,-30.)/2./PI); // 20.5 punt
   Eelectrons[10]=log10(pow(10,-30.)/2./PI); // 21.0 punt
   Eelectrons[11]=log10(pow(10,-30.)/2./PI); // 21.5 punt
   */


   // electron component of Figure 4 of ES&S
   // astro-ph/0101216
   Eelectrons[0] = -17.2; // 16.
   Eelectrons[1] = -17.35; // 16.5
   Eelectrons[2] = -17.2; // 17.
   Eelectrons[3] = -17.1; // 17.5
   Eelectrons[4] = -17.2; // 18.
   Eelectrons[5] = -17.5; // 18.5
   Eelectrons[6] = -18.0; // 19
   Eelectrons[7] = -18.5; // 19.5
   Eelectrons[8] = -19.4; // 20.
   Eelectrons[9] = -30.; // 20.5 punt
   Eelectrons[10] = -30.; // 21.0 punt
   Eelectrons[11] = -30.; // 21.5 punt
   // upper curve in Figure 9 of ES&S
   // astro-ph/0101216
   Emuons[0] = -16.85; //16.
   Emuons[1] = -16.4; //16.5
   Emuons[2] = -16.05; //17.
   Emuons[3] = -16.; // 17.5
   Emuons[4] = -16.15; // 18.
   Emuons[5] = -16.5; // 18.5
   Emuons[6] = -17.1; // 19
   Emuons[7] = -17.7; // 19.5
   Emuons[8] = -18.65; // 20.
   Emuons[9] = -19.75; // 20.5 punt
   Emuons[10] = -30.; // 21.0 punt
   Emuons[11] = -30.; // 21.5 punt


   for (int i = 0; i < 12; i++) {
      EdNdEdAdt[i] = pow(10., Eelectrons[i]) + pow(10., Emuons[i]);
   }

   maxflux = Max(EdNdEdAdt, 12);
   for (int i = 0; i < 12; i++) {
      EdNdEdAdt[i] = EdNdEdAdt[i] / maxflux;
   } //for

}

void GetGZK7()
{
   double Emuons[7];
   // E^2 dN/dE/dA/dt for neutrinos that are produced as muon neutrinos or muon antineutrinos.
   //KD I think it's EdN/dE/dA/dt
   double Eelectrons[7];// E^2 dN/dE/dA/dt for neutrinos that are produced as electron neutrinos or muon antineutrinos.
   //KD I think it's EdN/dE/dA/dt

   for (int i = 0; i < 7; i++) {
      energy[i] = 17. + ((double)i) / 2.;
      Emuons[i] = -30.;
      Eelectrons[i] = -30.;
   }

   // electron component of Figure 4 of ES&S
   // astro-ph/0101216
   Eelectrons[0] = -17.2; // 17.
   Eelectrons[1] = -17.1; // 17.5
   Eelectrons[2] = -17.2; // 18.
   Eelectrons[3] = -17.5; // 18.5
   Eelectrons[4] = -18.0; // 19
   Eelectrons[5] = -18.5; // 19.5
   Eelectrons[6] = -19.4; // 20.

   // upper curve in Figure 9 of ES&S
   // astro-ph/0101216
   Emuons[0] = -16.05; //17.
   Emuons[1] = -16.; // 17.5
   Emuons[2] = -16.15; // 18.
   Emuons[3] = -16.5; // 18.5
   Emuons[4] = -17.1; // 19
   Emuons[5] = -17.7; // 19.5
   Emuons[6] = -18.65; // 20.

   for (int i = 0; i < 7; i++) {
      EdNdEdAdt[i] = pow(10., Eelectrons[i]) + pow(10., Emuons[i]);
   }
   maxflux = Max(EdNdEdAdt, 7);
   for (int i = 0; i < 7; i++) {
      EdNdEdAdt[i] = EdNdEdAdt[i] / maxflux;
   } //for

}

double  PickEnergy(double* energy, double* EdNdEdAdt)
{

   double thisenergy = 16.;
   double thisflux = 2.;
   double max = 1.;
   int energybin = 0;
   while (thisflux > max) {
      // pick an energy
      thisenergy = Rand3.Rndm() * (Max(energy, 12) - Min(energy, 12));
      energybin = (int)(thisenergy / 0.5);
      max = EdNdEdAdt[energybin];
      thisflux = Rand3.Rndm();
   } //while
   EXPONENT = 16 + thisenergy;

   return pow(10., thisenergy + Min(energy, 12));
} //Get energy spectrum


void GetE2_fan()//from 10^17 to 10^20
{
   double E[] = { -14.7, -15.2, -15.7, -16.2, -16.7, -17.2, -17.7}; //from 17 to 20

   for (int i = 0; i < 7; i++) {
      energy[i] = 17. + ((double)i) / 2.;
   }
   for (int i = 0; i < 7; i++) {
      EdNdEdAdt[i] = pow(10., E[i]);
   }
   maxflux = Max(EdNdEdAdt, 7);
   for (int i = 0; i < 7; i++) {
      EdNdEdAdt[i] = EdNdEdAdt[i] / maxflux;
   }

}

void GetGZK_fan()
{
   double E[] = { -15.8, -15.6, -15.775, -16.125, -16.625, -17.25, -18.075}; //from 17 to 20

   for (int i = 0; i < 7; i++) {
      energy[i] = 17. + ((double)i) / 2.;

   }
   for (int i = 0; i < 7; i++) {
      EdNdEdAdt[i] = pow(10., E[i]);
   }
   maxflux = Max(EdNdEdAdt, 7);
   for (int i = 0; i < 7; i++) {
      EdNdEdAdt[i] = EdNdEdAdt[i] / maxflux;
   }
}



double  PickEnergy_fan(double* energy, double* EdNdEdAdt)
{
   double thisenergy = 17.;
   double thisflux = 2.;
   double max = 1.;
   int energybin = 0;
   while (thisflux > max) {
      // pick an energy
      thisenergy = Rand3.Rndm() * (Max(energy, 7) - Min(energy, 7));
      energybin = (int)(thisenergy / 0.5);
      max = EdNdEdAdt[energybin];
      thisflux = Rand3.Rndm();
   } //while
   EXPONENT = 17 + thisenergy;

   return pow(10., thisenergy + Min(energy, 7));

} //Get energy spectrum

//some functions for reading in from xml using tinyxml2
int SetIntValueXML(tinyxml2::XMLNode * pRoot, int &IntParam, const char* ParamName){

  tinyxml2::XMLElement * InputParam = pRoot->FirstChildElement(ParamName);
  if (!InputParam) {printf("Could not find %s. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  tinyxml2::XMLError eResult = InputParam->QueryIntText(&IntParam);
  if (eResult != tinyxml2::XML_SUCCESS) { printf("Coult not read %s as type int, tinyXML Error: %i\n",ParamName, eResult); return eResult; }
  printf("%s= %i\n",ParamName,IntParam);

  return 0;
}

int SetBoolValueXML(tinyxml2::XMLNode * pRoot, int &BoolParam, const char* ParamName){

  tinyxml2::XMLElement * InputParam = pRoot->FirstChildElement(ParamName);
  if (!InputParam) {printf("Could not find %s. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  tinyxml2::XMLError eResult = InputParam->QueryIntText(&BoolParam);
  if (eResult != tinyxml2::XML_SUCCESS) { printf("Coult not read %s as type int, tinyXML Error: %i\n",ParamName, eResult); return eResult; }
  if (BoolParam != 0 && BoolParam != 1) { printf("%s should be boolean, so it must be 0 or 1",ParamName); return 1; }
  printf("%s= %i\n",ParamName,BoolParam);

  return 0;
}

int SetDoubleValueXML(tinyxml2::XMLNode * pRoot, double &DoubleParam, const char* ParamName){

  tinyxml2::XMLElement * InputParam = pRoot->FirstChildElement(ParamName);
  if (!InputParam) {printf("Could not find %s. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  tinyxml2::XMLError eResult = InputParam->QueryDoubleText(&DoubleParam);
  if (eResult != tinyxml2::XML_SUCCESS) { printf("Coult not read %s as type double, tinyXML Error: %i\n",ParamName, eResult); return eResult; }
  printf("%s= %f\n",ParamName,DoubleParam);

  return 0;
}

int SetTextValueXML(tinyxml2::XMLNode * pRoot, const char * &TextParam, const char* ParamName){

  tinyxml2::XMLElement * InputParam = pRoot->FirstChildElement(ParamName);
  if (!InputParam) {printf("Could not find %s. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  tinyxml2::XMLText * textNode = InputParam->FirstChild()->ToText();
  if (!textNode) {printf("Could not convert %s to text. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  TextParam = textNode->Value();
  printf("%s= %s\n",ParamName,TextParam);

  return 0;
}

int SetTextValueXML(tinyxml2::XMLNode * pRoot, std::string &TextParam, const char* ParamName){

  tinyxml2::XMLElement * InputParam = pRoot->FirstChildElement(ParamName);
  if (!InputParam) {printf("Could not find %s. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  tinyxml2::XMLText * textNode = InputParam->FirstChild()->ToText();
  if (!textNode) {printf("Could not convert %s to text. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  TextParam = textNode->Value();
  printf("%s= %s\n",ParamName,TextParam.c_str());

  return 0;
}

int SetElementXML(tinyxml2::XMLNode * pRoot, tinyxml2::XMLElement * &Element, const char* ParamName){

  Element = pRoot->FirstChildElement(ParamName);
  if (!Element) {printf("Could not find %s. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  //  printf("Set pointer for %s= %s\n",ParamName,Element);

  return 0;
}

int SetElementXML(tinyxml2::XMLElement * Element, tinyxml2::XMLElement * &SubElement, const char* ParamName){
  SubElement = Element->FirstChildElement(ParamName);
  if (!SubElement) {printf("Could not find %s. tintXML Parsing Error: %i\n",ParamName,tinyxml2::XML_ERROR_PARSING_ELEMENT); return tinyxml2::XML_ERROR_PARSING_ELEMENT;}
  //  printf("Set pointer for %s= %s\n",ParamName,SubElement);

  return 0;
}

int ReadEventList(string inputFN, int &NNU, vector<EvtParams> &EventList){
  //This follows closely the  structure of the tinyxml2 tutorial found at shilohjames.wordpress.com/2014/04/27/tinyxml2-tutorial/
  //Sets NNU and EventList from xml file
  const char * infn = inputFN.c_str();
  cout<<"Reading Event List from "<<infn<<"..."<<endl;
  tinyxml2::XMLDocument xmlDoc;
  tinyxml2::XMLError eResult = xmlDoc.LoadFile(infn);
  if (eResult != tinyxml2::XML_SUCCESS) { printf("Error loading Event List file %s, tinyXML Error: %i\n",infn, eResult); return eResult; }

  tinyxml2::XMLNode * pRoot = xmlDoc.FirstChild();
  if (!pRoot) { printf("tinyXML Error reading input file: %i\n", tinyxml2::XML_ERROR_FILE_READ_ERROR); return tinyxml2::XML_ERROR_FILE_READ_ERROR; }

  int ErrorStat = 0;
  const char* ErrMesg = "Error in ReadEventList\n";

  //Now we just measure the size of VectAntennas, so remove a human derp point
  //ErrorStat = SetIntValueXML(pRoot, NAntPerStn, "N_Ant_perST");
  //if (ErrorStat) {printf(ErrMesg); return 1;}

  tinyxml2::XMLElement * pElement;
  ErrorStat = SetElementXML(pRoot,pElement,"List");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  tinyxml2::XMLElement * pListElement;
  ErrorStat = SetElementXML(pElement,pListElement,"Event");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  while (pListElement){
    EvtParams iEvent;

    ErrorStat = SetIntValueXML(pListElement, iEvent.EvID, "EvID");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetDoubleValueXML(pListElement, iEvent.Exponent, "Exponent");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetDoubleValueXML(pListElement, iEvent.Inelasticity, "Inelasticity");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetDoubleValueXML(pListElement, iEvent.X, "X");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetDoubleValueXML(pListElement, iEvent.Y, "Y");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetDoubleValueXML(pListElement, iEvent.Z, "Z");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetDoubleValueXML(pListElement, iEvent.ThetaDir, "ThetaDir");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetDoubleValueXML(pListElement, iEvent.PhiDir, "PhiDir");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetTextValueXML(pListElement, iEvent.Flavor, "Flavor");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetTextValueXML(pListElement, iEvent.Current, "Current");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    EventList.push_back(iEvent);
    pListElement = pListElement->NextSiblingElement("Event");
  }

  //Measure the size of VectAntennas to set NAntPerStn
  NNU = EventList.size();
  printf("NNU = %i\n",NNU);

  return 0;
}

int ReadStnGeo(string inputFN, int &NAntPerStn, vector<AntennaPlacement> &VectAntennas){
  //This follows closely the  structure of the tinyxml2 tutorial found at shilohjames.wordpress.com/2014/04/27/tinyxml2-tutorial/
  //Sets NAntPerStn and VectAntennas from xml file

  const char * infn = inputFN.c_str();

  cout<<"Reading Station Geometry from "<<infn<<"..."<<endl;

  tinyxml2::XMLDocument xmlDoc;
  tinyxml2::XMLError eResult = xmlDoc.LoadFile(infn);
  if (eResult != tinyxml2::XML_SUCCESS) { printf("Error loading Stn Geometry file %s, tinyXML Error: %i\n",infn, eResult); return eResult; }

  tinyxml2::XMLNode * pRoot = xmlDoc.FirstChild();
  if (!pRoot) { printf("tinyXML Error reading input file: %i\n", tinyxml2::XML_ERROR_FILE_READ_ERROR); return tinyxml2::XML_ERROR_FILE_READ_ERROR; }

  int ErrorStat = 0;
  const char* ErrMesg = "Error in ReadStnGeo\n";

  //Now we just measure the size of VectAntennas, so remove a human derp point
  //ErrorStat = SetIntValueXML(pRoot, NAntPerStn, "N_Ant_perST");
  //if (ErrorStat) {printf(ErrMesg); return 1;}

  tinyxml2::XMLElement * pElement;
  ErrorStat = SetElementXML(pRoot,pElement,"Antennas");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  tinyxml2::XMLElement * pListElement;
  ErrorStat = SetElementXML(pElement,pListElement,"Antenna");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  while (pListElement){
    AntennaPlacement iAntenna;

    ErrorStat = SetIntValueXML(pListElement, iAntenna.Type, "AntType");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    tinyxml2::XMLElement * pSubElement;
    ErrorStat = SetElementXML(pListElement,pSubElement,"position");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    tinyxml2::XMLElement * pValue;
    ErrorStat = SetElementXML(pSubElement,pValue,"value");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    for (int i =0; i<3; i++){
      eResult = pValue->QueryDoubleText(&iAntenna.position[i]);
      if (eResult != tinyxml2::XML_SUCCESS) { printf("Error reading position value. tinyXML Error: %i\n", eResult); return eResult; }
      pValue = pValue->NextSiblingElement("value");
    }

    ErrorStat = SetElementXML(pListElement,pSubElement,"n_boresight");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetElementXML(pSubElement,pValue,"value");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    for (int i =0; i<3; i++){
      eResult = pValue->QueryDoubleText(&iAntenna.n_boresight[i]);
      if (eResult != tinyxml2::XML_SUCCESS) { printf("Error reading n_boresight value. tinyXML Error: %i\n", eResult); return eResult; }
      pValue = pValue->NextSiblingElement("value");
    }

    ErrorStat = SetElementXML(pListElement,pSubElement,"n_epol");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    ErrorStat = SetElementXML(pSubElement,pValue,"value");
    if (ErrorStat) {printf(ErrMesg); return 1;}

    for (int i =0; i<3; i++){
      eResult = pValue->QueryDoubleText(&iAntenna.n_epol[i]);
      if (eResult != tinyxml2::XML_SUCCESS) { printf("Error reading n_epol value. tinyXML Error: %i\n", eResult); return eResult; }
      pValue = pValue->NextSiblingElement("value");
    }

    VectAntennas.push_back(iAntenna);
    pListElement = pListElement->NextSiblingElement("Antenna");
  }
  //Measure the size of VectAntennas to set NAntPerStn
  NAntPerStn = VectAntennas.size();
  printf("NAntPerStn = %i\n",NAntPerStn);

  return 0;
}

int ReadInputXML(const char * infn){
  //This follows closely the  structure of the tinyxml2 tutorial found at shilohjames.wordpress.com/2014/04/27/tinyxml2-tutorial/
  //Since this is ShelfMC, all of these parameters should be declaired at global scope ('cause that's the way we do things 'round here, apparently)

  cout<<"Reading Input Parameters from "<<infn<<"..."<<endl;

  tinyxml2::XMLDocument xmlDoc;
  tinyxml2::XMLError eResult = xmlDoc.LoadFile(infn);
  if (eResult != tinyxml2::XML_SUCCESS) { printf("Error loading input file, tinyXML Error: %i\n", eResult); return eResult; }

  tinyxml2::XMLNode * pRoot = xmlDoc.FirstChild();
  if (!pRoot) { printf("tinyXML Error reading input file: %i\n", tinyxml2::XML_ERROR_FILE_READ_ERROR); return tinyxml2::XML_ERROR_FILE_READ_ERROR; }

  int ErrorStat = 0;
  const char* ErrMesg = "Error in ReadInputXML\n";

  //Boolean parameters can go here (They are actually still ints, but who's watching?)
  ErrorStat = SetBoolValueXML(pRoot, USELIST, "USELIST");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, SPECTRUM, "SPECTRUM");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, WIDESPECTRUM, "WIDESPECTRUM");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, GZK, "GZK");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, FANFLUX, "FANFLUX");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, SIGNAL_FLUCT, "SIGNAL_FLUCT");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, TAUREGENERATION, "TAUREGENERATION");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, SHADOWING, "SHADOWING");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, H_PROP, "H_PROP");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, DEPTH_DEPENDENT_N, "DEPTH_DEPENDENT_N");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, FIRN, "FIRN");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, CONST_ATTENLENGTH, "CONST_ATTENLENGTH");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, HEXAGONAL, "HEXAGONAL");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  //The next few parameters are obsolete and should be removed
  ErrorStat = SetBoolValueXML(pRoot, SCATTER, "SCATTER");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetBoolValueXML(pRoot, DIPOLE, "DIPOLE");
  if (ErrorStat) {printf(ErrMesg); return 1;}


  //Int parameters can go here
  ErrorStat = SetIntValueXML(pRoot, NNU, "NNU");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetIntValueXML(pRoot, seed, "seed");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetIntValueXML(pRoot, NROWS, "NROWS");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetIntValueXML(pRoot, NCOLS, "NCOLS");
  if (ErrorStat) {printf(ErrMesg); return 1;}


  //N_Ant_perST is now set in the Stn Geometry File
  //  ErrorStat = SetIntValueXML(pRoot, N_Ant_perST, "N_Ant_perST");
  //  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetIntValueXML(pRoot, N_Ant_Trigger, "N_Ant_Trigger");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  //Double parameters can go here
  ErrorStat = SetDoubleValueXML(pRoot, EXPONENT, "EXPONENT");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, TNOISE, "TNOISE");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, ICETHICK, "ICETHICK");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, FIRNDEPTH, "FIRNDEPTH");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, NFIRN, "NFIRN");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, C_INDEX, "C_INDEX");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, ATTEN_UP, "ATTEN_UP");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, ATTEN_DOWN, "ATTEN_DOWN");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, ATTEN_H, "ATTEN_H");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, ATTEN_FACTOR, "ATTEN_FACTOR");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, REFLECT_RATE, "REFLECT_RATE");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, ATGap, "ATGap");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, ST4_R, "ST4_R");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, NSIGMA, "NSIGMA");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, gainv, "gainv");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, FREQ_LOW, "FREQ_LOW");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, FREQ_HIGH, "FREQ_HIGH");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, MAX_DISTANCE, "MAX_DISTANCE");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  //The next few parameters are obsolete and should be removed
   ErrorStat = SetDoubleValueXML(pRoot, SCATTER_WIDTH, "SCATTER_WIDTH");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  ErrorStat = SetDoubleValueXML(pRoot, Z, "Z");
  if (ErrorStat) {printf(ErrMesg); return 1;}

  //Text parameters (const char *) can go here
  ErrorStat = SetTextValueXML(pRoot, EventListFN, "EventListFN");
  if (ErrorStat) {printf("Error in ReadInputXML"); return 1;}

  ErrorStat = SetTextValueXML(pRoot, StnGeoFN, "StnGeoFN");
  if (ErrorStat) {printf("Error in ReadInputXML"); return 1;}

  return ErrorStat;
}
