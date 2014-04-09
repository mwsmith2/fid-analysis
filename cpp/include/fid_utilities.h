#ifndef __FID_UTILITIES_H__
#define __FID_UTILITIES_H__

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <vector>
#include <random>
using std::vector;
using std::string;

//--- root includes ---------------------------------------------------------//
#include "TCanvas.h"
#include "TGraph.h"

//--- project includes ------------------------------------------------------//
#include "fid_class.h"
#include "fid.h"

typedef vector<double> vec;

namespace fid{

  void DrawFID(const vec &wf, const vec &tm, 
    const string title, const string filename);

  void DrawFID(fid::FID &my_fid, const string title, const string filename);

  // Declare methods
  void AddWhiteNoise(vec &wf, double s2n=100.0);
  void ConstructTimeVector(int ntimes, double t0, double dt, vec &tm);
  void ConstructLinearGradient(int npoints, vec &grad);
  void ConstructQuadraticGradient(int npoints, vec &grad);
//  void ConstructGradientFID(vec &grad, vec &wf);

} // fid

#endif