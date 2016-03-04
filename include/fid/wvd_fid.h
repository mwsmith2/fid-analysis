#ifndef LIBFID_FID_WVD_FID_H_
#define LIBFID_FID_WVD_FID_H_

/*============================================================================*\

author: Ronaldo Ortez
email: supron00@gmail.com

notes:

  Class used to prepare normal FIDs for the calculation of the WVD, and additional functionality for the manipulation and plotting of WVD results.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include <complex>
#include <cmath>

//--- other includes --------------------------------------------------------//
#include <boost/algorithm/string.hpp>
#include <armadillo>
#include "TGraph.h"
#include "TF1.h"

//--- project includes ------------------------------------------------------//
#include "fid/base_fid.h"
#include "fid/params.h"
#include "fid/math.h"

namespace fid {

class WvdFid : public BaseFid {

 public:
  
  // ctors
  WvdFid(const std::string& fid_file, const bool upsample = false, const int window =0);
  WvdFid(const char* fid_file, const bool upsample = false, const int window = 0);
  WvdFid(const std::vector<double>& wf, const bool upsample = false, const int window = 0);
  WvdFid(const std::vector<double>& wf, const std::vector<double>& tm,
         const bool upsample = false, const int window = 0);
 
  // Member Functions
  double CalcFreq() {};// BaseFid virtual function inheritance override.
  // @Todo: Zero time frequency extraction
  // @Todo: WVD matrix extraction 
  // @Todo: Plotting functions
  // @Todo: Diagnostics

 private:
  
  int window_ ;
  bool upsample_;
  
  std::vector<double> wf_cen_;
  std::vector<double> wf_cen_up_;
  std::vector<cdouble> wf_analytic_;

  std::vector<cdouble> acf_;
  std::vector<double> wvd_f_;

  std::vector<double> HilbertRot(const std::vector<double>& v, const cdouble i = cdouble(0.0, 0.0));
  std::vector<cdouble> AutoCorr(int idx = 0);

  void WvdInit();
  void AnalyticFid();
  void WvdFreqExt();  
  void WvdCenter();

  void InitHook() {};//BaseFid virtual function inheritance override.
  
  
}; // WvdFid
 
} // fid

#endif
