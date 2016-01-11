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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace fid {

// Enumerate the different methods of frequency extraction
enum Method { ZC, CN, AN, LZ, EX, PH, SN,
              ZEROCOUNT, 
              CENTROID, 
              ANALYTICAL, 
              LORENTZIAN, 
              EXPONENTIAL, 
              PHASE,
              SINUSOID };

// constants
const double kTau = 2 * M_PI;

// general fid analysis parameters
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
  extern double snr_thresh;  // max_amp_ / noise_
  extern double len_thresh;  // fraction of signal

  extern Method freq_method;

} // ::params

// simulation parameters
namespace sim {

  extern int seed;   // the seed for the random generator
  extern double dt_integration; // step size of time for integration
  extern double snr;   // default signal-to-noise ratio
  extern double amplitude; // the generic amplitude of FID oscillations
  extern double baseline;  // the baseline offset for the FID waveform

  extern double start_time; // start time for the FID
  extern double delta_time; // time spacing for the FID
  extern int num_samples;    // the number of FID samples

  extern double freq_ref;    // reference frequency used to mix down
  extern double freq_larmor; // Larmor frequency to be simulated
  extern double freq_cut_ratio; // lowpass freq as fraction of freq_larmor;
  extern double mixdown_phi; // arbitrary phase for mixing freq
  extern std::vector<double> spin_0; // the initial spin vector

  extern double gamma_1;  // relaxation time gamma_1
  extern double gamma_2;  // relaxation time gamma_2
  extern double gamma_g;  // gyromagnetic ratio of proton

  extern double omega_r;  // strength of rf-pulse
  extern double t_pulse;  // length of nmr pulse

} // ::sim

namespace grad {

  extern std::string root_file;
  extern std::string fid_branch;
  extern double min;
  extern double max;
  extern int poln_order;
  extern std::vector<double> poln_coefs;
}

// header implementation of load_params
void load_params(std::string conf_file);

} // ::fid

#endif
