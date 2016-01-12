#ifndef LIBFID_INCLUDE_FID_FID_FACTORY_H_
#define LIBFID_INCLUDE_FID_FID_FACTORY_H_

/*===========================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu
file:  fid_sim.h

notes:

  The class, FidFactory simulates FIDs using the Bloch equations and 
  numerical integration.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <complex>

// //--- other includes --------------------------------------------------------//
#include <boost/numeric/odeint.hpp>
#include <armadillo>
#include "TFile.h"
#include "TTree.h"

// //--- project includes ------------------------------------------------------//
#include "fid/params.h"
#include "fid/math.h"

//--- namespaces ------------------------------------------------------------//

namespace fid {

class FidFactory
{
 public:

  // ctors
  FidFactory();

  // dtors
  ~FidFactory();

  // member methods
  void IdealFid(std::vector<double>& wf, std::vector<double>& tm);
  void SimulateFid(std::vector<double>& wf, std::vector<double>& tm);
  void GradientFid(const std::vector<double>& grad, std::vector<double>& wf);

  // Functions to set simulation parameters.
  void SetSeed(const int seed) { seed_ = seed; };
  void SetDtIntegration(const double step) { dt_integration_ = step; };

  void SetSNR(const double snr) { snr_ = snr; };
  void SetAmplitude(const double amp) { amplitude_ = amp; };
  void SetBaseline(const double baseline) { baseline_ = baseline; };
  void SetStartTime(const double start_time) { start_time_ = start_time; };
  void SetDeltaTime(const double delta_time) { delta_time_ = delta_time; };
  void SetNumSamples(int n) { num_samples_ = n; };

  void SetFreqLarmor(const double f) { freq_larmor_ = f; };
  void SetFreqRef(const double f) { freq_ref_ = f; };
  void SetFreqCutRatio(const double ratio) { freq_cut_ratio_ = ratio; };
  void SetMixdownPhi(const double phi) { mixdown_phi_ = phi; };
  void SetInitialSpin(const std::vector<double>& s) { spin_0_ = s; };

  void SetGamma1(const double g) { gamma_1_ = g; };
  void SetGamma2(const double g) { gamma_2_ = g; };
  void SetGammaG(const double g) { gamma_g_ = g; };
  void SetOmegaR(const double o) { omega_r_ = o; };
  void SetTPulse(const double t) { t_pulse_ = t; };

  void SetWithNoise(bool whith_noise) { with_noise_ = whith_noise; };
  void SetDiscretize(bool discretize) { discretize_ = discretize; };

  void PrintDiagnosticInfo();

 private:

  double ti_;
  double tf_;
  double dt_;
  int sim_to_fid_;
  int sim_length_;
  int printer_idx_;

  std::vector<double> s_;
  std::vector<double> spin_vec_;
  std::vector<double> time_vec_;
  std::vector<double> spin_sum_;
  std::vector<double> cos_cache_;
  std::vector<double> gradient_;

  // Simulation parameters (loaded from fid::sim namespace).
  int seed_;
  double dt_integration_;
  double snr_;
  double amplitude_;
  double baseline_;

  double start_time_;
  double delta_time_;
  int num_samples_;

  double freq_ref_;
  double freq_larmor_;
  double freq_cut_ratio_;
  double mixdown_phi_;
  std::vector<double> spin_0_;

  double gamma_1_;
  double gamma_2_;
  double gamma_g_;

  double omega_r_;
  double t_pulse_;

  bool discretize_;
  bool with_noise_;
  
  // Low pass filter to extract mixed down signal.
  std::vector<double> LowPassFilter(std::vector<double>& s);

  // Function which returns time dependent Bfield.
  std::vector<double> Bfield(const double& t);

  // The time evolution equation for the fields.
  void Bloch(std::vector<double> const &s, std::vector<double> &dsdt, double t);

  // The integration monitor function
  void Printer(std::vector<double> const &s , double t);

  int num_sim_fids_;
  int zero_idx_;
  double d_grad_;

  TFile *pf_fid_;
  TTree *pt_fid_;
  std::vector<Double_t> wf_;

  int GetTreeIndex(double grad_strength);

  void LoadParams();
}; // FidFactory

} // ::fid

#endif
