#ifndef _PCTCUT_H_
#define _PCTCUT_H_
// Class for selecting events to pass on to the image reconstruction
#include "TkrHits.h"
#include "pCT_Tracking.h"
#include "pCTconfig.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"

using namespace std;
class TkrHits;
class pCT_Tracking;
class pCTcut {
  int event_counter;
  int n1track; // Various counters for summarizing the number of events killed by cuts
  int nLT8hits;
  int nGoodXtra;
  int nBack;
  int nHitReject;
  int mxXlyr;
  int mxTotXlyr;
  int mxXhits;
  int Thread;
  int minTkrs;
  int maxTkrs;
  int pileUp; 

 public:
  static inline pCTcut* GetInstance() {return theCuts;} 
  int nKeep;
  // constructor
  pCTcut(); // Class constructor called prior to the event loop
  // classes
  bool cutHitSlope(int, int, float); 
  bool cutTrackIsocenterIntercept(float);
  void AddToPileUp();
  bool cutEvt(pCT_Tracking&, TkrHits&);
  bool EnrgCut(float [5], float, float, float, float );
  
//  bool dEEFilter(float, float, float*, float*);
  bool dEEFilter(float, float, TGraphErrors*);

  void CalculatedEEFilterParameters(TH2D* dEEhist, float dEElow[3], float dEEhigh[3], int stage);
  void summary();

  float mxSlope[2][2]; // Cut on the slope of the front tracker vector, separately for V and T
  float deltaMx;                    // Cut on how far the two vectors miss each other at u=0, in mm
 
 private: 
  static pCTcut *theCuts; 
  pCTconfig* theConfig;
};
#endif
