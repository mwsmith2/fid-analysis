#ifndef FID_ANALYSIS_INCLUDE_FID_UTILITIES_H_
#define FID_ANALYSIS_INCLUDE_FID_UTILITIES_H_

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
#include "fid_params.h"
#include "fid_class.h"
#include "fid.h"

namespace fid{

  void DrawFID(const vec &wf, const vec &tm, 
    const string title, const string filename);

  void DrawFID(fid::FID &my_fid, const string title, const string filename);

  // Declare methods
  void AddWhiteNoise(vec &wf, double snr=100.0);
  void ConstructTimeVector(int ntimes, double t0, double dt, vec &tm);
  void ConstructLinearGradient(int npoints, vec &grad);
  void ConstructQuadraticGradient(int npoints, vec &grad);
//  void ConstructGradientFID(vec &grad, vec &wf);

  template<typename T>
  inline vector<T> ConstructSweepRange(const vector<T>& range_vec){
    vector<T> res;

    for (T it = range_vec[0]; it <= range_vec[1]; it += range_vec[2]){
      res.push_back(it);
    }

    return res;
  }

} // fid

#endif