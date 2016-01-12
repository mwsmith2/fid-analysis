#ifndef LIBFID_INCLUDE_FID_PARAMS_EXTDEF_H_
#define LIBFID_INCLUDE_FID_PARAMS_EXTDEF_H_

#include "fid/params.h"

namespace fid {

std::string logdir("/var/log/fid/");

// general fid analysis parameters
namespace params {

  int fft_peak_width = 20;
  int zc_width = 100;
  int edge_ignore = 50;
  double start_thresh = 0.37;
  double max_phase_jump = 4.71;
  double low_pass_freq = 2000.0;
  double centroid_thresh = 0.01;
  double hyst_thresh = 0.7;
  double snr_thresh = 10.0;  // max_amp_ / noise_
  double len_thresh = 0.025;  // fraction of signal
  Method freq_method = PH;

} // ::params

// simulation parameters
namespace sim {

  int seed = 0;
  double dt_integration = 2.0e-5;
  double snr = 90000.0;   
  double amplitude = 2000.0;
  double baseline = 21000.0;

  int num_samples = 10000;
  double start_time = -1.0; 
  double delta_time = 0.001; 

  double freq_ref = 950.0;    
  double freq_larmor = 997.0; 
  double freq_cut_ratio = 0.1;
  double mixdown_phi = 0.0; 
  std::vector<double> spin_0 = {0.0, 0.0, 1.0};         

  double gamma_1 = 0.05;  
  double gamma_2 = 0.05;  
  double gamma_g = 1.0;  

  double omega_r = 50.0;  
  double t_pulse = 0.005; 

  bool with_noise = true;
  bool discretize = false;

} // ::sim

namespace grad {

  std::string root_file = "~/.fid/sim_fids.root";
  std::string fid_branch = "fid";
  double min = -2000;
  double max = 2000;
  int poln_order = 2;
  std::vector<double> poln_coefs = {0.0, 0.0, 1.0, 0.0};
}

// header implementation of load_params
void load_params(std::string conf_file)
{
  // using directives
  using boost::property_tree::ptree;
  using namespace sim;
  using namespace params;

  // Declare and load JSON params.
  ptree pt;
  read_json(conf_file, pt);
  pt = pt.get_child("fid");

  // general fid parameters
  num_samples = pt.get<int>("num_samples", num_samples);
  start_time = pt.get<double>("start_time", start_time);
  delta_time = pt.get<double>("delta_time", delta_time);  

  // analysis parameters
  fft_peak_width = pt.get<int>("params.fft_peak_width", fft_peak_width);
  edge_ignore = pt.get<int>("params.edge_ignore", edge_ignore);
  zc_width = pt.get<int>("params.zc_width", zc_width);
  start_thresh = pt.get<double>("params.start_thresh", start_thresh);
  max_phase_jump = pt.get<double>("params.max_phase_jump", max_phase_jump);
  low_pass_freq = pt.get<double>("params.low_pass_freq", low_pass_freq);
  centroid_thresh = pt.get<double>("params.centroid_thresh", centroid_thresh);
  hyst_thresh = pt.get<double>("params.hyst_thresh", hyst_thresh);

  // sim parameters
  seed = pt.get<int>("sim.seed", seed);
  dt_integration = pt.get<double>("sim.dt_integration", dt_integration);
  snr  = pt.get<double>("sim.snr", snr);

  gamma_1 = pt.get<double>("sim.gamma_1", gamma_1);
  gamma_2 = pt.get<double>("sim.gamma_2", gamma_2);
  gamma_g = pt.get<double>("sim.gamma_g", gamma_g);

  freq_larmor = pt.get<double>("sim.freq_larmor", freq_larmor);
  freq_ref = pt.get<double>("sim.freq_ref", freq_ref);
  freq_cut_ratio = pt.get<double>("sim.freq_cut_ratio", freq_cut_ratio);
  mixdown_phi = pt.get<double>("sim.mixdown_phi", mixdown_phi);

  omega_r = pt.get<double>("sim.omega_r", omega_r);
  t_pulse = pt.get<double>("sim.t_pulse", t_pulse);
  with_noise = pt.get<bool>("sim.with_noise", with_noise);
  discretize = pt.get<bool>("sim.discretize", discretize);

  // gradient fid file parameters
  grad::root_file = pt.get<std::string>("grad.root_file", grad::root_file);
  grad::fid_branch = pt.get<std::string>("grad.fid_branch", grad::fid_branch);
  grad::min = pt.get<double>("grad.min", grad::min);
  grad::max = pt.get<double>("grad.max", grad::max);
  grad::poln_order = pt.get<double>("grad.poln_order", grad::poln_order);

  try {
    std::vector<double> tmp;
    for (auto &v : pt.get_child("grad.poln_coefs")){
      tmp.push_back(v.second.get_value<double>());
    }
    grad::poln_coefs = tmp;
  }
  catch (boost::property_tree::ptree_bad_path){};

} // load_params

} // ::fid

#endif
