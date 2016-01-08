#include "fid/utils.h"
#include "fid/params_extdef.h"

namespace fid {

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
  fit_width = pt.get<int>("params.fit_width", fit_width);
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


// Create a vector of linear spaced values.
std::vector<double> linspace(double x0, double xf, int n)
{
  // If n isn't set assuming integer spacing.
  if (n == 0) {
    n = xf - x0 + 0.5;
  }

  double dx = (xf - x0) / n;
  std::vector<double> vec(0.0, n);

  for (int i = 0; i < n; ++i){
    vec[i] = dx * i + x0;
  }

  return vec;
}


// Create a discrete quadratic gradient with a normalized unit amplitude.
std::vector<double> normalized_gradient(int npoints, int poln)
{
  std::vector<double> grad;

  // First get the spacing right.
  for (int i = 0; i < npoints; i++){
    grad.push_back(pow(i, poln));
  }

  // Subtract off the average.
  double avg = std::accumulate(grad.begin(), grad.end(), 0.0) / grad.size();

  for (uint i = 0; i < grad.size(); i++){
    grad[i] -= avg;
  }

  // Normalize by largest value.
  double max = *std::max_element(grad.begin(), grad.end());

  for (uint i = 0; i < grad.size(); i++){
    grad[i] /= max;
  }

  return grad;
}

} // ::fid
