#ifndef pedGainCalib_h
#define pedGainCalib_h

// Class to encapsulate all of Vladimir Bashkirov pedestal and gain on-the-fly calibration
// R.P. Johnson   5/22/2016

#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

#include "TFile.h"
#include "TVcorrection.h"
#include "TH1D.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "pCTraw.h"
#include "pCTconfig.h"
#include "TFitResult.h"
using namespace std;

class TFitResult;

class pedGainCalib {
 public:
 pedGainCalib(TFile* root, int pedMin[5], float oldPed[5], float t1, float t2, float t3, float t4);
 
  ~pedGainCalib();  
  //Variables
  float Ped[5];
  float stdPed[5];
  float FivepercentPed[5];
  float Peak[5];
  float stdGain[5];
  TH1D* hPed[5];
  TH2D* hPed2D[5];
  long long unsigned int startTime; 
  bool first = true;
  TH1D* hPed_In[5];
  TH1D* hPed_Out[5];
  TH1D* hTotFil[5];
  TH1D* hTotUnFil[5];
  TH1D* hEnrg[5];
  TH1D* hEnrgTot;
  TProfile2D* hProfT;
  TH1D* hTedet;
  TFile* RootFile;
  float emtrng1, emtrng2, emtrng3, emtrng4;
  string inFileName_s;
  float GainFac[5]; // Gain correction factor derived here for each stage

  // functions
  void ClearHist();
  void FillPeds(pCTraw &rawEvt);
  void GetPeds();
  void FillADC(int ADC[5]);
  void ResetHist();
  void WriteHist();
  void FillGains(float V[4], float T[4], float Ene[5], int phSum[5]);
  void GetGains(TVcorrection *TVcorr);
 private:
  pCTconfig* theConfig;
  
}; // end of class pedGainCalib
#endif
