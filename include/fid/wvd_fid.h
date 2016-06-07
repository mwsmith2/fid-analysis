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
#include "TMultiGraph.h"
#include "TF1.h"
#include "TH1D.h"

//--- project includes ------------------------------------------------------//
#include "fid/base_fid.h"
#include "fid/params.h"
#include "fid/math.h"

namespace fid {

class WvdFid : public BaseFid {

 public:
  
  // ctors
    WvdFid(const std::string& fid_file, const int fid_stop = 0, const int window = 80, const bool upsample = true);
    WvdFid(const char* fid_file, const int fid_stop = 0, const int window = 80, const bool upsample = true);
    WvdFid(const std::vector<double>& wf, const int fid_stop = 0, const int window = 80, const bool upsample = true);
    WvdFid(const std::vector<double>& wf, const std::vector<double>& tm, const int fid_stop =0,
         const int window =80, const bool upsample = true);
  
  // accessors
  const std::vector<double>& wvd_f() const {return wvd_f_; };
  const std::vector<cdouble>& AnalyticWf() const {return  wf_analytic_;};
  const std::vector<cdouble>& acf() const {return acf_;};

  // Member Functions
  double CalcFreq() {};// BaseFid virtual function inheritance override.
  void WvdFreqExt();
  double WvdZeroFreq(const int fid_stop = 1500);
  // @Todo: Zero time frequency extraction
  // @Todo: WVD matrix extraction 
  // @Todo: Plotting functions
  // @Todo: Diagnostics
 
  
  // Plotting Functions
  void SaveWvd(std::string filename, std::string title= "");
  void MutltiPlot(const int pads = 0);

 private:
  
  int window_ ;
  int fid_stop_ ;
  double dt_;
  double wvd_f0_;
  bool upsample_ ;
  
  std::vector<double> wf_cen_;
  std::vector<double> wf_cen_up_;
  std::vector<cdouble> wf_analytic_;

  std::vector<cdouble> acf_;
  std::vector<double> wvd_f_;

  std::vector<cdouble> Filter(const std::vector<cdouble> &v, const int freq_cut=0);
  std::vector<double> HilbertRot(const std::vector<double>& v, const cdouble i = cdouble(0.0, 0.0), const int freq_cut =0);
  std::vector<cdouble> AutoCorr(const std::vector<cdouble> &wf_temp, int idx = 0);

  void WvdInit();
  void AnalyticFid(); 
  void WvdCenter();

  void InitHook() {};//BaseFid virtual function inheritance override.
  
  
}; // WvdFid
 
} // fid

#endif
