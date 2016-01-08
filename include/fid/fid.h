#ifndef FID_ANALYSIS_FID_FID_H_
#define FID_ANALYSIS_FID_FID_H_

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
#include "TCanvas.h"
#include "TF1.h"

//--- project includes ------------------------------------------------------//
#include "fid/params.h"
#include "fid/math.h"
#include "fid/utils.h"

namespace fid {

class FID {

 public:
  
  // ctors
  FID(const std::string& fid_file);
  FID(const char* fid_file);
  FID(const std::vector<double>& wf, const std::vector<double>& tm);
  FID(const std::vector<double>& wf);

  // Simplified frequency extraction
  double GetFreq();
  double GetFreq(const std::string& method_name);
  double GetFreq(const char* method_name);
  double GetFreq(const Method m);
  double GetFreqError();
  
  // specific frequency extraction methods
  double CalcFreq();
  
  // diagnostic function
  void PrintDiagnosticInfo();
  void PrintDiagnosticInfo(std::iostream out);

  // utility functions
  void SaveData(std::string filename);
  void SavePlot(std::string filename, std::string title="");
  void SaveTimeFit(std::string filename, std::string title="");
  void SaveFreqFit(std::string filename, std::string title="");
  void SaveTimeRes(std::string filename, std::string title="");
  void SaveFreqRes(std::string filename, std::string title="");

  // accessors
  const bool isgood() const { return health_ > 0.0; };
  const std::vector<double>& wf() const { return wf_; };
  const std::vector<double>& tm() const { return tm_; };
  const std::vector<double>& res() const { return res_; };
  const std::vector<double>& power() const { return power_; };
  const std::vector<double>& fftfreq() const { return fftfreq_; };
  const std::vector<double>& phase() const { return phase_ ;}
  const std::vector<double>& env() const { return env_; };
  const double& freq() const {  return freq_;  }
  const double& freq_err() const { return freq_err_; };
  const double& chi2() const { return chi2_; };
  const double &health() const { return health_; };
  const double fid_time() const { return tm_[f_wf_ - 1] - tm_[i_wf_]; };
  const double snr() const { return pow(max_amp_ / noise_, 2); };
  const TGraph& gr_time_series() const { return gr_time_series_; };
  const TGraph& gr_freq_series() const { return gr_freq_series_; };
  const TF1&    f_fit() const { return f_fit_; };
  const unsigned int& i_wf() { return i_wf_; };
  const unsigned int& f_wf() { return f_wf_; };
  const unsigned int& i_fft() { return i_fft_; };
  const unsigned int& f_fft() { return f_fft_; };

  // More specific use functions.
  double CalcZeroCountFreq();
  double CalcCentroidFreq();
  double CalcAnalyticalFreq();
  double CalcLorentzianFreq();
  double CalcSoftLorentzianFreq();
  double CalcExponentialFreq();
  double CalcPhaseFreq(int poln=1);
  double CalcPhaseDerivFreq(int poln=1);
  double CalcSinusoidFreq();
  
 private:
  
  // Private Member Variables
  unsigned int i_wf_; // start and stop of relevant data
  unsigned int f_wf_;
  unsigned int i_fft_;
  unsigned int f_fft_;
  double health_;
  double noise_;
  double max_amp_;
  double mean_;
  double chi2_; // Store the most recent chi2
  double freq_;
  double freq_err_;
  Method freq_method_;

  std::vector<double> guess_;
  TF1 f_fit_;
  TGraph gr_time_series_;
  TGraph gr_freq_series_;

  // For general plotting.
  TCanvas c1_;
  TGraph gr_;
  
  // bigger data arrays
  std::vector<double> wf_;
  std::vector<double> tm_;
  std::vector<double> res_;
  std::vector<double> power_;
  std::vector<double> env_;
  std::vector<double> phase_;
  std::vector<double> fftfreq_;
  std::vector<double> temp_; // for random transformations

  // Private Member Functions  
  // init function to be called after wf_ and tm_ are set.
  void Init();
  void SaveGraph(std::string filename, std::string title);
  void LoadTextData(std::string filename); 

  // internal utility functions
  void CalcNoise();
  void CalcMaxAmp();			
  void CenterFid();
  void FindFidRange();
  void CalcPowerEnvAndPhase();
  void CalcFftFreq();
  void GuessFitParams();
  void FreqFit(TF1& func);
  Method ParseMethod(const std::string& m);
  
}; // FID
 
} // fid

#endif
