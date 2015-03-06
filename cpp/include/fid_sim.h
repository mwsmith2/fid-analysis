#ifndef FID_ANALYSIS_INCLUDE_FID_SIM_H_
#define FID_ANALYSIS_INCLUDE_FID_SIM_H_

/*===========================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu
file:  fid_sim.h

notes:

  The class, FidFactory simulates FIDs using the Bloch equations and 
  numerical integration.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <algorithm>
#include <random>
#include <cmath>
#include <functional>

//--- other includes --------------------------------------------------------//
#include <boost/numeric/odeint.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <armadillo>
#include "TFile.h"
#include "TTree.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_math.h"

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
  void IdealFid(std::vector<double>& wf, 
                std::vector<double>& tm, 
                bool withnoise=false);

  void SimulateFid(std::vector<double>& wf, 
                   std::vector<double>& tm, 
                   bool withnoise=false);

  void GradientFid(const std::vector<double>& gradient, 
                   std::vector<double>& wf, 
                   bool withnoise=false);

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
}; // FidFactory

} // ::fid

#endif