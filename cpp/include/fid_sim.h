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

using namespace std::chrono;
using namespace boost::numeric::odeint;
using namespace boost::property_tree;

namespace fid {

class FidFactory
{
public:

  // ctors
  FidFactory();

  // dtors
  ~FidFactory();

  // member methods
  void IdealFid(vec& wf, vec& tm);
  void SimulateFid(vec& wf, vec& tm);
  void GradientFid(const vec& gradient, vec& wf);
  void PrintDiagnosticInfo();

 private:

  double ti_;
  double tf_;
  double dt_;
  int sim_to_fid_;
  int sim_length_;
  int printer_idx_;

  vec s_;
  vec spin_vec_;
  vec time_vec_;
  vec spin_sum_;
  vec cos_cache_;
  vec gradient_;
  
  // Low pass filter to extract mixed down signal.
  vec LowPassFilter(vec& s);

  // Function which returns time dependent Bfield.
  vec Bfield(const double& t);

  // The time evolution equation for the fields.
  void Bloch(vec const &s, vec &dsdt, double t);

  // The integration monitor function
  void Printer(vec const &s , double t);

  int num_sim_fids_;
  int zero_idx_;
  double d_grad_;

  TFile *pf_fid_;
  TTree *pt_fid_;
  vector<Double_t> wf_;

  int GetTreeIndex(double grad_strength);
}; // FidFactory

} // ::fid

#endif