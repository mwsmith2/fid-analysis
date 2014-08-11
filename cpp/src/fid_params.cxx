#include "fid_params.h"

namespace fid {

// per configuration statistics
int num_fids;
int len_fids;

// time vector variables
double i_time;
double d_time;

// default parameters
double s_freq;
double s_phase;
double s_grad;
double s_snr;
double s_tau;

// sweep parameters
namespace sweep {

  // freqeuency sweep
  bool freq_sweep;
  vec  freq_range;

  // phase sweep
  bool phase_sweep;
  vec  phase_range;

  // gradient sweep
  bool grad_sweep;
  vec  grad_range;

  // signal-to-noise ratio sweep
  bool snr_sweep;
  vec  snr_range;

} // ::sweep

// simulation parameters
namespace sim {

  int seed;
  vec spin_0;

  double dt_integration;

  double omega_r;
  double t_pulse;
  double gamma_g;
  double gamma_1;
  double gamma_2;
  double freq_ref;
  double freq_larmor;

} // ::sim

namespace grad {
  string root_file;
  string fid_branch;
  double min;
  double max;
  int poln_order;
  vec poln_coefs;
}

// general fid analysis params
namespace params {

  int fit_width;
  int edge_ignore;
  int zc_width;
  double start_thresh;
  double zc_alpha;
  double max_phase_jump;
  double low_pass_freq;
  double centroid_thresh;
  double hyst_thresh;

} // ::params

void load_params(int argc, char **argv)
{
  // using directives
  using boost::property_tree::ptree;
  using namespace sweep;
  using namespace sim;
  using namespace params;

  // load defaults first
  string config_file("runtime/.default_fid_params.json");
  ptree pt;
  read_json(config_file, pt);
  pt = pt.get_child("fid");

  // general fid parameters
  num_fids = pt.get<int>("num_fids");
  len_fids = pt.get<int>("len_fids");
  i_time = pt.get<double>("i_time");
  d_time = pt.get<double>("d_time");  

  s_freq = pt.get<double>("s_freq");
  s_phase  = pt.get<double>("s_phase");
  s_snr  = pt.get<double>("s_snr");
  s_grad = pt.get<double>("s_grad");
  s_tau = pt.get<double>("s_tau");

  // sweep parameters
  freq_sweep = pt.get<bool>("sweep.freq_sweep");
  BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.freq_range")){
    freq_range.push_back(v.second.get_value<double>());
  }

  phase_sweep = pt.get<bool>("sweep.phase_sweep");
  BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.phase_range")){
    phase_range.push_back(v.second.get_value<double>());
  }

  grad_sweep = pt.get<bool>("sweep.grad_sweep");
  BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.grad_range")){
    grad_range.push_back(v.second.get_value<double>());
  }

  snr_sweep = pt.get<bool>("sweep.snr_sweep");
  BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.snr_range")){
    snr_range.push_back(v.second.get_value<double>());
  }

  // sim parameters
  seed = pt.get<int>("sim.seed");
  BOOST_FOREACH(ptree::value_type &v, pt.get_child("sim.spin_0")){
    spin_0.push_back(v.second.get_value<double>());
  }

  dt_integration = pt.get<double>("sim.dt_integration");

  omega_r = pt.get<double>("sim.omega_r");
  t_pulse = pt.get<double>("sim.t_pulse");
  gamma_g = pt.get<double>("sim.gamma_g");
  gamma_1 = pt.get<double>("sim.gamma_1");
  gamma_2 = pt.get<double>("sim.gamma_2");
  freq_ref = pt.get<double>("sim.freq_ref");
  freq_larmor = pt.get<double>("sim.freq_larmor");

  // gradient fid file parameters
  grad::root_file = pt.get<string>("grad.root_file");
  grad::fid_branch = pt.get<string>("grad.fid_branch");
  grad::min = pt.get<double>("grad.min");
  grad::max = pt.get<double>("grad.max");
  grad::poln_order = pt.get<double>("grad.poln_order");

  BOOST_FOREACH(ptree::value_type &v, pt.get_child("grad.poln_coefs")){
    grad::poln_coefs.push_back(v.second.get_value<double>());
  }

  // analysis parameters
  fit_width = pt.get<int>("params.fit_width");
  edge_ignore = pt.get<int>("params.edge_ignore");
  zc_width = pt.get<int>("params.zc_width");
  start_thresh = pt.get<double>("params.start_thresh");
  zc_alpha = pt.get<double>("params.zc_alpha");
  max_phase_jump = pt.get<double>("params.max_phase_jump");
  low_pass_freq = pt.get<double>("params.low_pass_freq");
  centroid_thresh = pt.get<double>("params.centroid_thresh");
  hyst_thresh = pt.get<double>("params.hyst_thresh");

  // If the user provided a different config file, load it instead
  if (argc < 2) {
    config_file = string("./runtime/fid_params.json");
  } else {
    config_file = string(argv[1]);
  }
  read_json(config_file, pt);
  pt = pt.get_child("fid");

  // general fid parameters
  num_fids = pt.get<int>("num_fids", num_fids);
  len_fids = pt.get<int>("len_fids", len_fids);
  i_time = pt.get<double>("i_time", i_time);
  d_time = pt.get<double>("d_time", d_time);  

  s_freq = pt.get<double>("s_freq", s_freq);
  s_phase  = pt.get<double>("s_phase", s_phase);
  s_snr  = pt.get<double>("s_snr", s_snr);
  s_grad = pt.get<double>("s_grad", s_grad);
  s_tau = pt.get<double>("s_tau", s_tau);

  // sweep parameters
  freq_sweep = pt.get<bool>("sweep.freq_sweep", freq_sweep);
  try {
    vec tmp;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.freq_range")){
      tmp.push_back(v.second.get_value<double>());
    }
    freq_range = tmp;
  }
  catch (boost::property_tree::ptree_bad_path){};

  phase_sweep = pt.get<bool>("sweep.phase_sweep");
  try {
    vec tmp;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.phase_range")){
      tmp.push_back(v.second.get_value<double>());
    }
    phase_range = tmp;
  }
  catch (boost::property_tree::ptree_bad_path){};

  grad_sweep = pt.get<bool>("sweep.grad_sweep");
  try {
    vec tmp;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.grad_range")){
      tmp.push_back(v.second.get_value<double>());
    }
    grad_range = tmp;
  }
  catch (boost::property_tree::ptree_bad_path){};

  snr_sweep = pt.get<bool>("sweep.snr_sweep");
  try {
    vec tmp;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("sweep.snr_range")){
      tmp.push_back(v.second.get_value<double>());
    }
    snr_range = tmp;
  }
  catch (boost::property_tree::ptree_bad_path){};

  // sim parameters
  seed = pt.get<int>("sim.seed", seed);

  try {
    vec tmp;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("sim.spin_0")){
      tmp.push_back(v.second.get_value<double>());
    }
    spin_0 = tmp;
  }
  catch (boost::property_tree::ptree_bad_path){};

  dt_integration = pt.get<double>("sim.dt_integration", dt_integration);

  omega_r = pt.get<double>("sim.omega_r", omega_r);
  t_pulse = pt.get<double>("sim.t_pulse", t_pulse);
  gamma_g = pt.get<double>("sim.gamma_g", gamma_g);
  gamma_1 = pt.get<double>("sim.gamma_1", gamma_1);
  gamma_2 = pt.get<double>("sim.gamma_2", gamma_2);
  freq_ref = pt.get<double>("sim.freq_ref", freq_ref);
  freq_larmor = pt.get<double>("sim.freq_larmor", freq_larmor);

  // gradient fid file parameters
  grad::root_file = pt.get<string>("grad.root_file", grad::root_file);
  grad::fid_branch = pt.get<string>("grad.fid_branch", grad::fid_branch);
  grad::min = pt.get<double>("grad.min", grad::min);
  grad::max = pt.get<double>("grad.max", grad::max);
  grad::poln_order = pt.get<double>("grad.poln_order", grad::poln_order);

  try {
    vec tmp;
    BOOST_FOREACH(ptree::value_type &v, pt.get_child("grad.poln_coefs")){
      tmp.push_back(v.second.get_value<double>());
    }
    grad::poln_coefs = tmp;
  }
  catch (boost::property_tree::ptree_bad_path){};

  // analysis parameters
  fit_width = pt.get<int>("params.fit_width", fit_width);
  edge_ignore = pt.get<int>("params.edge_ignore", edge_ignore);
  zc_width = pt.get<int>("params.zc_width", zc_width);
  start_thresh = pt.get<double>("params.start_thresh", start_thresh);
  zc_alpha = pt.get<double>("params.zc_alpha", zc_alpha);
  max_phase_jump = pt.get<double>("params.max_phase_jump", max_phase_jump);
  low_pass_freq = pt.get<double>("params.low_pass_freq", low_pass_freq);
  centroid_thresh = pt.get<double>("params.centroid_thresh", centroid_thresh);
  hyst_thresh = pt.get<double>("params.hyst_thresh", hyst_thresh);

} // load_params

} // ::fid
