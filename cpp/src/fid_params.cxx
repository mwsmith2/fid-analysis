#include "fid_params.h"


namespace fid {

// per configuration statistics
int num_fids;
int len_fids;

// time vector variables
double i_time;
double d_time;

// default parameters
double freq;
double phi;
double grad;
double snr;

// sweep parameters
namespace sweep {

  // freqeuency sweep
  double i_freq;
  double f_freq;
  double d_freq;

  // phase sweep
  int num_phi;
  double i_phi;
  double f_phi;
  double d_phi;

  // gradient sweep
  double i_grad;
  double f_grad;
  double d_grad;

  // signal-to-noise sweep
  double i_snr;
  double f_snr;
  double d_snr;
  
} // ::sweep

// simulation parameters
namespace sim {

  int num_points;
  int reduction;
  int num_steps;

  double d_bfield;
  double i_time;
  double f_time;
  double dt;

  double omega_r;
  double t_pulse;
  double gamma_g;
  double gamma_1;
  double gamma_2;

} // ::sim

// general fid analysis params
namespace params {

  int fit_width;
  int zc_width;
  int ph_edge_ignore;
  double start_thresh;
  double zc_alpha;
  double ph_max_jump;

} // ::params

void load_params(int argc, char **argv)
{
  // load defaults first
  string config_file("runtime/.default_fid_params.json");
  read_json(config_file, pt);

  // general fid parameters
  num_fids = pt.get<int>("num_fids");
  len_fids = pt.get<int>("len_fids");
  i_time = pt.get<double>("i_time");
  d_time = pt.get<double>("d_time");  

  freq = pt.get<double>("freq");
  phi  = pt.get<double>("phi");
  snr  = pt.get<double>("snr");
  grad = pt.get<double>("grad");

  // sweep parameters
  sweep::i_freq = get.get<double>("sweep.i_freq");
  sweep::f_freq = get.get<double>("sweep.f_freq");
  sweep::d_freq = get.get<double>("sweep.d_freq");

  sweep::num_phi = get.get<int>("sweep.num_phi");
  sweep::i_phi = get.get<double>("sweep.i_phi");
  sweep::f_phi = get.get<double>("sweep.f_phi");

  sweep::i_grad = get.get<double>("sweep.i_grad");
  sweep::f_grad = get.get<double>("sweep.f_grad");
  sweep::d_grad = get.get<double>("sweep.d_grad");

  sweep::i_snr = get.get<double>("sweep.i_snr");
  sweep::f_snr = get.get<double>("sweep.f_snr");
  sweep::d_snr = get.get<double>("sweep.d_snr");

  // sim parameters
  sim::num_points = get.get<int>("sim.num_points");
  sim::num_steps = get.get<int>("sim.num_steps");
  sim::reduction = get.get<int>("sim.reduction");

  sim::d_bfield = get.get<double>("sim.d_bfield");
  sim::dt_integration = get.get<double>("sim.dt_integration");

  sim::omega_r = get.get<double>("sim.omega_r");
  sim::t_pulse = get.get<double>("sim.t_pulse");
  sim::gamma_g = get.get<double>("sim.gamma_g");
  sim::gamma_1 = get.get<double>("sim.gamma_1");
  sim::gamma_2 = get.get<double>("sim.gamma_2");
  sim::freq_ref = get.get<double>("sim.freq_ref");
  sim::freq_larmor = get.get<double>("sim.freq_larmor");

  // analysis parameters
  params::fit_width = get.get<int>("params.fit_width");
  params::start_thresh = get.get<int>("params.start_thresh");
  params::zc_width = get.get<int>("params.zc_width");
  params::zc_alpha = get.get<int>("params.zc_alpha");
  params::ph_edge_ignore = get.get<int>("params.ph_edge_ignore");
  params::ph_max_jump = get.get<int>("params.ph_max_jump");

  // If the user provided a different config file, load it instead
  if (argc < 2) return;
  config_file = string(argv[1]);

  // general fid parameters
  num_fids = pt.get_optional<int>("num_fids");
  len_fids = pt.get_optional<int>("len_fids");
  i_time = pt.get_optional<double>("i_time");
  d_time = pt.get_optional<double>("d_time");  

  freq = pt.get_optional<double>("freq");
  phi  = pt.get_optional<double>("phi");
  snr  = pt.get_optional<double>("snr");
  grad = pt.get_optional<double>("grad");

  // sweep parameters
  sweep::i_freq = pt.get_optional<double>("sweep.i_freq");
  sweep::f_freq = pt.get_optional<double>("sweep.f_freq");
  sweep::d_freq = pt.get_optional<double>("sweep.d_freq");

  sweep::num_phi = pt.get_optional<int>("sweep.num_phi");
  sweep::i_phi = pt.get_optional<double>("sweep.i_phi");
  sweep::f_phi = pt.get_optional<double>("sweep.f_phi");

  sweep::i_grad = pt.get_optional<double>("sweep.i_grad");
  sweep::f_grad = pt.get_optional<double>("sweep.f_grad");
  sweep::d_grad = pt.get_optional<double>("sweep.d_grad");

  sweep::i_snr = pt.get_optional<double>("sweep.i_snr");
  sweep::f_snr = pt.get_optional<double>("sweep.f_snr");
  sweep::d_snr = pt.get_optional<double>("sweep.d_snr");

  // sim parameters
  sim::num_points = pt.get_optional<int>("sim.num_points");
  sim::num_steps = pt.get_optional<int>("sim.num_steps");
  sim::reduction = pt.get_optional<int>("sim.reduction");

  sim::d_bfield = pt.get_optional<double>("sim.d_bfield");
  sim::dt_integration = pt.get_optional<double>("sim.dt_integration");

  sim::omega_r = pt.get_optional<double>("sim.omega_r");
  sim::t_pulse = pt.get_optional<double>("sim.t_pulse");
  sim::gamma_g = pt.get_optional<double>("sim.gamma_g");
  sim::gamma_1 = pt.get_optional<double>("sim.gamma_1");
  sim::gamma_2 = pt.get_optional<double>("sim.gamma_2");
  sim::freq_ref = pt.get_optional<double>("sim.freq_ref");
  sim::freq_larmor = pt.get_optional<double>("sim.freq_larmor");

  // analysis parameters
  params::fit_width = pt.get_optional<int>("params.fit_width");
  params::start_thresh = pt.get_optional<int>("params.start_thresh");
  params::zc_width = pt.get_optional<int>("params.zc_width");
  params::zc_alpha = pt.get_optional<int>("params.zc_alpha");
  params::ph_edge_ignore = pt.get_optional<int>("params.ph_edge_ignore");
  params::ph_max_jump = pt.get_optional<int>("params.ph_max_jump");

} // load_params

} // ::fid
