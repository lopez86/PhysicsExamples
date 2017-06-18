/*
  ROOT script to show how we can point to a position in the sky based on
  gravitational wave detector timing information.

  Calculates the chi^2 for a signal being at a particular location in the sky
  given the time difference.
  (c) Jeremy Lopez, 2016
*/
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <cmath>
//Coordinates in standard spherical coords (i.e. z-axis is north pole)
double hanford[3] = {0.7578,2.0837,1e-6};
double livingston[3] = {1.0386,1.5839,1e-6};
double virgo[3] = {0.8093,0.1833,1e-6};

double earthRadius = 6371e3;//m
double speedOfLight = 3e8;//m/s

TH2*
ChiSquare(double t, double unc=1e-6)
{
  double uncertainty=sqrt(2)*unc;
/*
  if (useTimeUncertainty){
    uncertainty = hypot(hanford[2],livingston[2]);
  }
*/
  double hf[3] = {sin(hanford[0])*cos(hanford[1]),sin(hanford[0])*sin(hanford[1]),cos(hanford[0])};
  double l[3] = {sin(livingston[0])*cos(livingston[1]),sin(livingston[0])*sin(livingston[1]),cos(livingston[0])};


  TH2F* chisq = new TH2F("chisq","chisq;#theta [rad];#phi [rad]",100,0,TMath::Pi(),100,0,2*TMath::Pi());
  for (int i =1 ; i <=100; i++)
  {
    double theta = chisq->GetXaxis()->GetBinCenter(i);
    for (int j = 1; j<=100; j++){
      double phi = chisq->GetYaxis()->GetBinCenter(j);
      double u[3] = {sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)};
      double hfCos = hf[0]*u[0]+hf[1]*u[1]+hf[2]*u[2];
      double lCos = l[0]*u[0]+l[1]*u[1]+l[2]*u[2];

      double dx = earthRadius * (hfCos-lCos);//effective spatial difference
      double dt = dx / speedOfLight;
      double chisqVal = (dt - t)*(dt-t) / (uncertainty*uncertainty);

      chisq->SetBinContent(i,j,exp(-0.5*chisqVal));
    }

  }  

  return chisq;

}
