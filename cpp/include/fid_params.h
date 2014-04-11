#ifndef FID_ANALYSIS_INCLUDE_FID_PARAMS_H_
#define FID_ANALYSIS_INCLUDE_FID_PARAMS_H_

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

//--- other includes --------------------------------------------------------//
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

namespace fid {

// using directives
using std::cout;
using std::vector;
using std::endl;
using std::string;
using boost::property_tree::ptree;

// typedefs
typedef vector<double> vec;

// constants
const double kTau = 2 * M_PI;

// per configuration statistics
extern int num_fids;
extern int len_fids;

// time vector variables
extern double i_time;
extern double d_time;

// default parameters
extern double s_freq;
extern double s_phase;
extern double s_grad;
extern double s_snr;

// sweep parameters
namespace sweep {

  // freqeuency sweep
  extern bool freq_sweep;
  extern vec  freq_range;

  // phase sweep
  extern bool phase_sweep;
  extern vec  phase_range;

  // gradient sweep
  extern bool grad_sweep;
  extern vec  grad_range;

  // signal-to-noise ratio sweep
  extern bool snr_sweep;
  extern vec  snr_range;

} // ::sweep

// simulation parameters
namespace sim {

  extern int num_points;
  extern int reduction;
  extern int num_steps;

  extern double d_bfield;
  extern double dt_integration;

  extern double omega_r;
  extern double t_pulse;
  extern double gamma_g;
  extern double gamma_1;
  extern double gamma_2;
  extern double freq_ref;
  extern double freq_larmor;

} // ::sim

// general fid analysis params
namespace params {

  extern int fit_width;
  extern int zc_width;
  extern int edge_ignore;
  extern double start_thresh;
  extern double zc_alpha;
  extern double max_phase_jump;
  extern double low_pass_freq;
  extern double centroid_thresh;
  extern double hyst_thresh;

} // ::params

void load_params(int argc, char **argv);

} // ::fid

#endif
