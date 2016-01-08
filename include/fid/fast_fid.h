#ifndef FID_ANALYSIS_FID_FAST_FID_H_
#define FID_ANALYSIS_FID_FAST_FID_H_

/*============================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes:

  Defines several frequency extraction and analysis methods for FIDs 
  as a class to encapsulate all the ideas.

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
#include "fid/params.h"
#include "fid/math.h"

namespace fid {

class FastFid {

 public:
  
  // ctors
  FastFid(const std::string& fid_file);
  FastFid(const char* fid_file);
  FastFid(const std::vector<double>& wf, const std::vector<double>& tm);
  FastFid(const std::vector<double>& wf);

  // Simplified frequency extraction
  double CalcFreq();
  double GetFreq();
  double GetFreqError();
  
  // diagnostic function
  void PrintDiagnosticInfo();
  void PrintDiagnosticInfo(std::iostream out);

  // accessors
  const bool isgood() const { return health_ > 0.0; };
  const std::vector<double>& wf() const { return wf_; };
  const std::vector<double>& tm() const { return tm_; };
  const double& freq() const {  return freq_;  }
  const double& freq_err() const { return freq_err_; };
  const ushort &health() const { return health_; };
  const double fid_time() const { return tm_[f_wf_] - tm_[i_wf_]; };
  const double snr() const { return pow(max_amp_ / noise_, 2); };
  const TGraph& gr_time_series() const { return gr_time_series_; };
  const unsigned int& i_wf() { return i_wf_; };
  const unsigned int& f_wf() { return f_wf_; };

 private:
  
  // Private Member Variables
  unsigned int i_wf_; // start and stop of relevant data
  unsigned int f_wf_;
  ushort health_; // percentage between 0 and 100.
  double noise_;
  double max_amp_;
  double mean_;
  double freq_;
  double freq_err_;
  Method freq_method_;

  std::vector<double> guess_;
  TGraph gr_time_series_;
  TGraph gr_freq_series_;
  
  // bigger data arrays
  std::vector<double> wf_;
  std::vector<double> tm_;
  std::vector<double> temp_; // for random transformations

  // Private Member Functions  
  // init function to be called after wf_ and tm_ are set.
  void Init();
  void LoadTextData(std::string filename);

  // internal utility functions
  void CalcNoise();
  void CalcMaxAmp();      
  void CenterFid();
  void FindFidRange();
  
}; // FastFid
 
} // fid

#endif
