#ifndef FID_ANALYSIS_INCLUDE_FID_PARAMS_H_
#define FID_ANALYSIS_INCLUDE_FID_PARAMS_H_

/*---------------------------------------------------------------------------*\

author: Matthias W. Smith
email: mwsmith2@uw.edu

about: This header file holds the projects parameter namespace.  These 
      variables are all either constants or initialized from json 
      configuration files.  This enables parameter tweaking without 
      recompilation.

\*---------------------------------------------------------------------------*/

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

// typedefs
typedef vector<double> vec;

// constants
const double kTau = 2 * M_PI;

// per configuration statistics
extern int num_fids; // the number of FIDs to analyze
extern int len_fids; // the number of samples in each FID

// time vector variables
extern double i_time; // start time for the FID
extern double d_time; // time spacing for the FID

// default parameters (s_ for single)
extern double s_freq;  // default frequency
extern double s_phase; // default phase
extern double s_grad;  // default gradient strength
extern double s_snr;   // default signal-to-noise ratio
extern double s_tau;   // default FID decay time

// sweep parameters
namespace sweep {

  // NOTE: ranges taken to be defined as [min, max, step]
  extern bool freq_sweep;  // sweep through frequencies
  extern vec  freq_range;  // define the desired frequency range
  extern bool phase_sweep; // sweep through phases
  extern vec  phase_range; // define the desired phase range
  extern bool grad_sweep;  // sweep through gradient strengths
  extern vec  grad_range;  // define the desired gradient strength range
  extern bool snr_sweep;   // sweep through signal-to-noise values
  extern vec  snr_range;   // define the signal-to-noise range

} // ::sweep

// simulation parameters
namespace sim {

  extern int seed;   // the seed for the random generator
  extern double dt_integration; // step size of time for integration

  extern vec spin_0;      // the initial spin vector
  extern double omega_r;  // strength of rf-pulse
  extern double t_pulse;  // length of nmr pulse
  extern double gamma_g;  // gyromagnetic ratio of proton
  extern double gamma_1;  // relaxation time gamma_1
  extern double gamma_2;  // relaxation time gamma_2
  extern double freq_ref; // reference frequency used to mix down
  extern double mixdown_phi; // arbitrary phase for mixing freq
  extern double freq_larmor; // Larmor frequency to be simulated

} // ::sim

namespace grad {

  extern string root_file;
  extern string fid_branch;
  extern double min;
  extern double max;
  extern int poln_order;
  extern vec poln_coefs;
}

// general fid analysis params
namespace params {

  extern int fit_width;  // fit width used by spectral peak fits
  extern int zc_width;   // size of window used to calculate FID noise
  extern int edge_ignore; // samples to ignore when doing phase fits
  extern double start_thresh; // threshold above noise to define start of FID
  extern double zc_alpha; // parameter used by exponential moving average
  extern double max_phase_jump; // maximum change allowed when unwrapping phase
  extern double low_pass_freq; // low pass frequency used by FID
  extern double centroid_thresh; // threshold of values included in centroid
  extern double hyst_thresh; // hysteris threshold used for zero counting

} // ::params

// run this function first thing in any module using this library
void load_params(int argc, char **argv);

} // ::fid

#endif
