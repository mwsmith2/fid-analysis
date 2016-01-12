#ifndef LIBFID_FID_BASE_FID_H_
#define LIBFID_FID_BASE_FID_H_

/*============================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes:

  Defines a base class for FID classes.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <string>
#include <vector>
#include <numeric>
#include <cmath>

//--- other includes --------------------------------------------------------//
#include "boost/filesystem.hpp"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"

//--- project includes ------------------------------------------------------//
#include "fid/params.h"
#include "fid/math.h"

namespace fid {

class BaseFid {

 public:

  // Default ctors/dtors.
  BaseFid() {};
  ~BaseFid() {};
  
  // Simplified frequency extraction
  virtual double CalcFreq() = 0;
  double GetFreq();
  double GetFreqError();
  
  // diagnostic function
  void PrintDiagnosticData(std::ostream& out=std::cout);
  void DumpDiagnosticData(std::string dirname=logdir, 
                          std::string filestub="");

  // accessors
  const std::vector<double>& wf() const { return wf_; };
  const std::vector<double>& tm() const { return tm_; };
  const double& freq() const {  return freq_;  }
  const double& freq_err() const { return freq_err_; };

  const double fid_time() const { return tm_[f_wf_] - tm_[i_wf_]; };
  const double snr() const { return pow(max_amp_ / noise_, 2); };
  const bool isgood() const { return health_ > 0.0; };
  const ushort &health() const { return health_; };

  const unsigned int& i_wf() { return i_wf_; };
  const unsigned int& f_wf() { return f_wf_; };
  const unsigned int& i_fft() { return i_fft_; };
  const unsigned int& f_fft() { return f_fft_; };

  const double& chi2() const { return chi2_; };
  const TF1&    f_fit() const { return f_fit_; };
  const TGraph& gr_time_series() const { return gr_time_series_; };
  const TGraph& gr_freq_series() const { return gr_freq_series_; };

  // Utility functions
  void SaveData(std::string filename);
  void SaveGraph(std::string filename, std::string title);
  void SavePlot(std::string filename, std::string title="");
  void SaveTimeFit(std::string filename, std::string title="");
  void SaveFreqFit(std::string filename, std::string title="");
  void SaveTimeRes(std::string filename, std::string title="");
  void SaveFreqRes(std::string filename, std::string title="");

 protected:
  
  // Private Member Variables
  unsigned int i_wf_; // start and stop of relevant data
  unsigned int f_wf_;
  unsigned int i_fft_;
  unsigned int f_fft_;

  // Waveform characteristics
  double mean_;
  double noise_;
  double max_amp_;

  double freq_;
  double freq_err_;
  double chi2_; // Store the most recent chi2
  ushort health_; // percentage between 0 and 100.

  // Load default (or user configured params)
  double edge_width_;
  double edge_ignore_;
  double start_amplitude_;
  double max_phase_jump_;
  double low_pass_freq_;
  double fft_peak_width_;
  double centroid_thresh_;
  double hyst_thresh_;
  double snr_thresh_;
  double len_thresh_;
  Method freq_method_;

  // For fits.
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
  std::vector<double> temp_; // for random transformations

  // Private Member Functions  
  // init function to be called after wf_ and tm_ are set.
  void Init();

  // internal utility functions
  void LoadTextData(std::string filename);
  void LoadParams();
  void CalcNoise();
  void CalcMaxAmp();      
  void CenterFid();
  void FindFidRange();
  void virtual InitHook() {};
  
}; // BaseFid
 
} // fid

#endif
