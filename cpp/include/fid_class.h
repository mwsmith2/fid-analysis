#ifndef FID_ANALYSIS_FID_CLASS_H_
#define FID_ANALYSIS_FID_CLASS_H_

/*===========================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes:

  This library consists of several frequency extraction and analysis 
  methods for FIDs as well as a class to encapsulate all the ideas.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>

//--- other includes --------------------------------------------------------//
#include <fftw3.h>
#include "TGraph.h"
#include "TF1.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_math.h"

namespace fid
{

class FID {

 public:
  
  // ctor
  FID(const vec& wf, const vec& tm);

  // Diagnostic Function
  void PrintDiagnosticInfo();
  void PrintDiagnosticInfo(std::iostream out);
  
  // frequency extraction methods
  double CalcZeroCountFreq();
  double CalcCentroidFreq();
  double CalcAnalyticalFreq();
  double CalcLorentzianFreq();
  double CalcSoftLorentzianFreq();
  double CalcExponentialFreq();
  double CalcPhaseFreq(int poln=1);
  double CalcPhaseDerivFreq(int poln=1);
  double CalcSinusoidFreq();
  
  // accessors
  const bool& isgood() const {return isgood_;};
  const vec& wf() const {return wf_;};
  const vec& tm() const {return tm_;};
  const vec& res() const {return res_;};
  const vec& power() const {return power_;};
  const vec& freq() const {return freq_;};
  const vec& phase() const {return phase_;}
  const vec& env() const {return env_;};
  const double& chi2() const {return chi2_;};
  const double& freq_err() const {return freq_err_;};
  const double snr() const {return max_amp_*max_amp_ / (noise_ * noise_);};
  const TGraph& gr_time_series() const {return gr_time_series_;};
  const TGraph& gr_freq_series() const {return gr_freq_series_;};
  const TF1&    f_fit() const {return f_fit_;};
  
 private:
  
  // Member Variables
  bool isgood_;
  uint i_wf_; // start and stop of relevant data
  uint f_wf_;
  uint i_fft_;
  uint f_fft_;
  double noise_;
  double max_amp_;
  double mean_;
  double chi2_; // Store the most recent chi2
  double freq_err;
  vec guess_;
  TF1 f_fit_;
  TGraph gr_time_series_;
  TGraph gr_freq_series_;
  
  // bigger data arrays
  vec wf_;
  vec tm_;
  vec res_;
  vec power_;
  vec env_;
  vec phase_;
  vec freq_;
  vec temp_; // for random transformations
  
  // internal utility functions
  void CalcNoise();
  void CalcMaxAmp();			
  void CenterFid();
  void FindFidRange();
  void CalcPowerEnvAndPhase();
  void CalcFftFreq();
  void GuessFitParams();
  void FreqFit(TF1& func);
  
}; // FID
 
} // fid

#endif
