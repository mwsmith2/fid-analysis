#include "fid_utils.h"
#include "fid_params_extdef.h"

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


// Fill the input waveform vector with an ideal FID.
void ideal_fid(std::vector<double>& wf, 
               std::vector<double>& tm, 
               double f, double phi, 
               double snr, 
               double tau, 
               double t0)
{
	wf.reserve(tm.size());
	wf.resize(0);

	// Define the waveform
	double temp;

	for (auto it = tm.begin(); it != tm.end(); it++){

		if (*it >= t0){
			temp = std::exp(-(*it - t0) / tau);
			temp *= std::sin((*it) * 2 * M_PI * f + phi);
			wf.push_back(temp);

		} else {
			wf.push_back(0.0);

		}
	} 

	// Add some noise
	static std::default_random_engine gen(0);
	std::normal_distribution<double> norm(0.0, 1.0 / snr);
	for (auto it = wf.begin(); it != wf.end(); it++){
		*it += norm(gen);
	}
}


// Create a linearly spaced vector over the interval.
void construct_time_vector(uint num_times, 
                           double t0, 
                           double dt, 
                           std::vector<double> &tm)
{
  if (tm.size() != num_times){
    tm.resize(num_times);
  }

  for (uint i = 0; i < num_times; i++){
    tm[i] = dt * i + t0;
  }
}


// Create a discrete quadratic gradient with a normalized unit amplitude.
void construct_quadratic_gradient(int npoints, std::vector<double> &grad)
{
  // First get the spacing right.
  for (int i = 0; i < npoints; i++){
    grad.push_back((double)i * i);
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
}


// Create a discrete linear gradient with a normalized unit amplitude.
void construct_linear_gradient(int npoints, std::vector<double> &grad)
{
  // First get the spacing right.
  for (int i = 0; i < npoints; i++){
    grad.push_back((double)i);
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
}


// Create a TCanvas and print the graph.
void draw_graph(TGraph gr, std::string fname, std::string title)
{
  // Set up the graph.
  std::string new_title(title);
  new_title.append("; time [ms]; amplitude [a.u.]");
  gr.SetTitle(new_title.c_str());

  // Get our own TCanvas
  TCanvas c1;

  // Draw the waveform
  gr.Draw();
  c1.Print(fname.c_str());
}


// Create a TGraph and feed it into the draw_graph() function.
void draw_graph(const std::vector<double> &wf, const std::vector<double> &tm, std::string fname, std::string title)
{
  // Create a graph
  TGraph gr(wf.size(), &tm[0], &wf[0]);

  draw_graph(gr, fname, title);
}


// Feed the waveform and time vectors into draw_graph.
void draw_fid(const FID &my_fid, std::string fname, std::string title)
{
  // Get the data vectors
  std::vector<double> wf = my_fid.wf();
  std::vector<double> tm = my_fid.tm();

  draw_graph(wf, tm, fname, title);
}


// Print the time series fit from an FID.
void draw_fid_time_fit(const FID &my_fid, std::string fname, std::string title)
{
  // Copy the graph
  TGraph gr = my_fid.gr_time_series();
  draw_graph(gr, fname, title);
}


// Print the time series residuals from an FID.
void draw_fid_time_res(const FID &my_fid, std::string fname, std::string title)
{
  // Copy the residuals
  std::vector<double> res = my_fid.res();

  // Copy the time series graph
  TGraph gr_fit = my_fid.gr_time_series();
  TGraph gr(res.size());

  // Set the points
  for (uint i = 0; i < res.size(); ++i){
    static double x, y;

    gr_fit.GetPoint(i, x, y);
    gr.SetPoint(i, x, res[i]); 
  }

  draw_graph(gr, fname, title);
}


// Print the freq series fit from an FID.
void draw_fid_freq_fit(const FID &my_fid, std::string fname, std::string title)
{
  // Copy the graph
  TGraph gr = my_fid.gr_freq_series();
  draw_graph(gr, fname, title);
}


// Print the freq series residuals from an FID.
void draw_fid_freq_res(const FID &my_fid, std::string fname, std::string title)
{
  // Copy the residuals
  std::vector<double> res = my_fid.res();

  // Copy the frequency series graph
  TGraph gr_fit = my_fid.gr_freq_series();
  TGraph gr(res.size());

  // Set the points
  for (uint i = 0; i < res.size(); ++i){
    static double x, y;

    gr_fit.GetPoint(i, x, y);
    gr.SetPoint(i, x, res[i]); 
  }

  draw_graph(gr, fname, title);
}


// Analyze an FID with all methods and print to output stream csv style.
void calc_freq_write_csv(FID& my_fid, std::ofstream& out)
{
  // Test all the frequency extraction methods and write the results
  out << my_fid.CalcZeroCountFreq() << ", " << 0.0 << ", ";
  out << my_fid.CalcCentroidFreq() << ", " << 0.0 << ", ";
  out << my_fid.CalcAnalyticalFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcLorentzianFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSoftLorentzianFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcExponentialFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSinusoidFreq() << ", " << my_fid.chi2() << std::endl;
}

// Analyze an FID with all methods and print to output stream csv style.
void calc_phase_freq_write_csv(FID& my_fid, std::ofstream& out)
{
  // Test all the frequency extraction methods using the phase
  out << my_fid.CalcPhaseFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseFreq(2) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseFreq(3) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseDerivFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseDerivFreq(2) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseDerivFreq(3) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSinusoidFreq() << ", " << my_fid.chi2() << std::endl;
}


// Read a FID text file into the waveform and time vectors.
void read_fid_file(std::string fname, 
                   std::vector<double> &wf, 
                   std::vector<double> &tm)
{
  // open the file first
  std::ifstream in(fname);

  // shrink vectors
  wf.resize(0);
  tm.resize(0);

  double wf_temp;
  double tm_temp;

  while (in.good())
  {
    in >> tm_temp >> wf_temp;
    tm.push_back(tm_temp);
    wf.push_back(wf_temp);
  }
}


// Write FID waveform/time vectors to a text file.
void write_fid_file(std::string fname,
                    const std::vector<double> &wf, 
                    const std::vector<double> &tm)
{
  // open the file first
  std::ofstream out(fname);

  for (int i = 0; i < tm.size(); ++i) {
    out << tm[i] << ", " << wf[i] << std::endl;
  }
}

} // ::fid
