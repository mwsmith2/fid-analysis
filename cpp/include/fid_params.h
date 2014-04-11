#ifndef FID_ANALYSIS_INCLUDE_FID_PARAMS_H_
#define FID_ANALYSIS_INCLUDE_FID_PARAMS_H_

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

namespace fid {

// using directives
using std::cout;
using std::vector;
using std::endl;

// typedefs
typedef vector<double> vec;

// constants
extern const double kTau = 2 * M_PI;

// sweep parameters
namespace sweep {

  // per configuration statistics
  extern int n_fids;

  // default parameters (s_ for single)
  extern double s_freq;
  extern double s_phi;
  extern double s_grad;
  extern double s_s2n;

  // time vector variables
  extern int fid_length;
  extern double i_time;
  extern double d_time;

  // freqeuency sweep
  extern double i_freq;
  extern double f_freq;
  extern double d_freq;

  // phase sweep
  extern int n_phi;
  extern double i_phi;
  extern double f_phi;
  extern double d_phi;

  // gradient sweep
  extern double i_grad;
  extern double f_grad;
  extern double d_grad;

  // signal-to-noise sweep
  extern double i_s2n;
  extern double f_s2n;
  extern double d_s2n;
  
} // ::params

// simulation parameters
namespace sim {

  extern int num_fids;
  extern int num_points;
  extern int reduction;
  extern int num_steps;

  extern double d_bfield;
  extern double ti;
  extern double tf;
  extern double dt;

  extern double omega_r;
  extern double t_pulse;
  extern double gamma_g;
  extern double gamma_1;
  extern double gamma_2;

} // ::sim

// general fid analysis params
extern int fit_width;
extern int zc_width;
extern int ph_edge_ignore;
extern double start_thresh;
extern double zc_alpha;
extern double ph_max_jump;

void load_params(int argc, char **argv);

} // ::fid
