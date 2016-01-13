#ifndef LIBFID_INCLUDE_FID_PARAMS_EXTDEF_H_
#define LIBFID_INCLUDE_FID_PARAMS_EXTDEF_H_

#include "fid/params.h"

namespace fid {

std::string logdir("/var/log/fid/");

// general fid analysis parameters
namespace params {

  double fft_peak_width = 20;
  double edge_width = 100;
  double edge_ignore = 50;
  double start_amplitude = 0.37;
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
  double integration_step = 2.0e-5;
  double signal_to_noise = 90000.0;   
  double amplitude = 2000.0;
  double baseline = 21000.0;

  int num_samples = 10000;
  double start_time = -1.0; 
  double sample_time = 0.001; 

  double mixdown_freq = 950.0;    
  double larmor_freq = 997.0; 
  double lowpass_ratio = 0.1;
  double mixdown_phi = 0.0; 
  std::vector<double> spin_0 = {0.0, 0.0, 1.0};         

  double gamma_1 = 0.05;  
  double gamma_2 = 0.05;  
  double gamma_g = 1.0;  

  double rf_omega = 50.0;  
  double rf_duration = 0.005; 

  bool addnoise = true;
  bool discrete = false;

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
  sample_time = pt.get<double>("sample_time", sample_time);  

  // analysis parameters
  fft_peak_width = pt.get<int>("params.fft_peak_width", fft_peak_width);
  edge_ignore = pt.get<int>("params.edge_ignore", edge_ignore);
  edge_width = pt.get<int>("params.edge_width", edge_width);
  start_amplitude = pt.get<double>("params.start_amplitude", start_amplitude);
  max_phase_jump = pt.get<double>("params.max_phase_jump", max_phase_jump);
  low_pass_freq = pt.get<double>("params.low_pass_freq", low_pass_freq);
  centroid_thresh = pt.get<double>("params.centroid_thresh", centroid_thresh);
  hyst_thresh = pt.get<double>("params.hyst_thresh", hyst_thresh);

  // sim parameters
  seed = pt.get<int>("sim.seed", seed);
  integration_step = pt.get<double>("sim.integration_step", integration_step);
  signal_to_noise  = pt.get<double>("sim.signal_to_noise", signal_to_noise);

  gamma_1 = pt.get<double>("sim.gamma_1", gamma_1);
  gamma_2 = pt.get<double>("sim.gamma_2", gamma_2);
  gamma_g = pt.get<double>("sim.gamma_g", gamma_g);

  larmor_freq = pt.get<double>("sim.larmor_freq", larmor_freq);
  mixdown_freq = pt.get<double>("sim.mixdown_freq", mixdown_freq);
  lowpass_ratio = pt.get<double>("sim.lowpass_ratio", lowpass_ratio);
  mixdown_phi = pt.get<double>("sim.mixdown_phi", mixdown_phi);

  rf_omega = pt.get<double>("sim.rf_omega", rf_omega);
  rf_duration = pt.get<double>("sim.rf_duration", rf_duration);
  addnoise = pt.get<bool>("sim.addnoise", addnoise);
  discrete = pt.get<bool>("sim.discrete", discrete);

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
