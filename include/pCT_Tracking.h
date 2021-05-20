#ifndef pCT_Tracking_h
#define pCT_Tracking_h

// Tracker pattern recognition routines
// R.P. Johnson  5/22/2016

#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>

#include "TkrHits.h"
#include "pCTgeo.h"
#include "pCTcut.h"
#include "pCTconfig.h"

struct Tkr2D {
  float X[4], U[4]; // the lateral (X) and longitudinal (U) coordinates of the track
  int Qh[4];         // hit quality: 0=interpolated, 1=beam constraint, 2=gap, 3=measured
  float Miss;       // distance in mm by which the two halves of the track miss each  other at u=0
  bool Good;         // true if the track candidate has not been rejected in favor of a different one
  int Q;             // track quality: 3= best, made from 4 measured points
                     //                2= made from 3 points and a T layer gap
                     //                1= made from 3 points and a beam origin constraint
                     //                0= made from 3 points and interpolation or extrapolation to
                     // the back layers
};

struct vct {
  float Y[2], U[2], I[2], intercept, slope;
};

class pCTcut;
class pCTconfig;
class pCT_Tracking {

  //Substitute Variables to aid readability

  float Y0_candidate,U0; // Y can mean either V or T here (lateral direction)
  float Y1_candidate,U1;
  float Y2_candidate,U2;
  float Y3_candidate,U3;
  float miss;

  inline float quadExtrap(float x1, float x2, float x3, float x4, float y1, float y2, float y3) {
    float m11 = x2 * x2 - x1 * x1;
    float m12 = x2 - x1;
    float m21 = x3 * x3 - x2 * x2;
    float m22 = x3 - x2;
    float v1 = y2 - y1;
    float v2 = y3 - y2;
    float det = m11 * m22 - m12 * m21;
    if (det == 0.) {
      std::cout << "pCT_Tracking::quadExtrap warning: a zero determinate was "
                   "encountered in making a quadratic exrapolation." << std::endl;
      return -9999.;
    }
    float a = (m22 * v1 - m12 * v2) / det;
    float b = (m11 * v2 - m21 * v1) / det;
    float c = y1 - (a * x1 + b) * x1;
    //        std::cout << std::endl;
    //        std::cout << y1 << " vs " << (a*x1 + b)*x1 + c << std::endl;
    //        std::cout << y2 << " vs " << (a*x2 + b)*x2 + c << std::endl;
    //        std::cout << y3 << " vs " << (a*x3 + b)*x3 + c << std::endl;
    //        std::cout << std::endl;
    return (a * x4 + b) * x4 + c;
  };

  // This is the pattern recognition, which works only in 2D, separately for the
  // V-U and T-U views.
  std::vector<Tkr2D> Tracking2D(int Idx, TkrHits &pCThits, pCTgeo* Geometry);

public:

  // To use these results for image reconstruction, require nTracks==1 and then
  // access the good tracks using itkV and itkT.
  int nTracks;                // Number of tracks found (max of V and T views if at least 1 in
                              // each, otherwise 0).
  std::vector<Tkr2D> VTracks; // List of all the tracks in the V view
  std::vector<Tkr2D> TTracks; // List of all the tracks in the T view
  int itkV;                   // Index of the first good track in the V view
  int itkT;                   // Index of the first good track in the T view

  inline float frontPredV(int tk, float u) { // Extrapolate the front V track
                                               // vector to the plane u
    float slope = (VTracks[tk].X[1] - VTracks[tk].X[0]) / (VTracks[tk].U[1] - VTracks[tk].U[0]);
    return VTracks[tk].X[1] + slope * (u - VTracks[tk].U[1]);
  }
  inline float frontPredT(int tk, float u) { // Extrapolate the front T track
                                               // vector to the plane u
    float slope = (TTracks[tk].X[1] - TTracks[tk].X[0]) / (TTracks[tk].U[1] - TTracks[tk].U[0]);
    return TTracks[tk].X[1] + slope * (u - TTracks[tk].U[1]);
  }
  inline float backPredV(int tk, float u) { // Extrapolate the back V track
                                              // vector to the plane u
    float slope = (VTracks[tk].X[3] - VTracks[tk].X[2]) / (VTracks[tk].U[3] - VTracks[tk].U[2]);
    return VTracks[tk].X[3] + slope * (u - VTracks[tk].U[3]);
  }
  inline float backPredT(int tk, float u) { // Extrapolate the back T track
                                              // vector to the plane u
    float slope = (TTracks[tk].X[3] - TTracks[tk].X[2]) / (TTracks[tk].U[3] - TTracks[tk].U[2]);
    return TTracks[tk].X[3] + slope * (u - TTracks[tk].U[3]);
  }

  pCT_Tracking(TkrHits &pCThits, pCTgeo* Geometry);

  // Method to print out the list of tracks
  void dumpTracks(int eventNumber);

private:
  pCTcut* theCuts; 
  pCTconfig* theConfig;

}; // End of the pCT_Tracking class
#endif
