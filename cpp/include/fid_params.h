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
extern double freq;
extern double phi;
extern double grad;
extern double snr;

// sweep parameters
namespace sweep {

  // freqeuency sweep
  extern double i_freq;
  extern double f_freq;
  extern double d_freq;

  // phase sweep
  extern int num_phi;
  extern double i_phi;
  extern double f_phi;
  extern double d_phi;

  // gradient sweep
  extern double i_grad;
  extern double f_grad;
  extern double d_grad;

  // signal-to-noise sweep
  extern double i_snr;
  extern double f_snr;
  extern double d_snr;
  
} // ::params

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
