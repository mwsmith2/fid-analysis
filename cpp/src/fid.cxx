#include "fid.h"

namespace fid
{
// Generate an ideal FID
void ideal_fid(vec& wf, vec& tm, double f, double phi, 
	double snr, double tau, double t0){

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

	return;
}

void construct_time_vector(uint num_times, double t0, double dt, vec &tm)
{
  if (tm.size() != num_times){
    tm.resize(num_times);
  }

  for (uint i = 0; i < num_times; i++){
    tm[i] = dt * i + t0;
  }
}


void construct_quadratic_gradient(int num_points, vec &grad)
{
  // construct a normalized, centered quadratic gradient

  // first get the spacing right
  for (int i = 0; i < num_points; i++){
    grad.push_back((double)i * i);
  }

  // subtract off the average
  double avg = std::accumulate(grad.begin(), grad.end(), 0.0) / grad.size();
  for (uint i = 0; i < grad.size(); i++){
    grad[i] -= avg;
  }

  // normalize by largest value
  double max = *std::max_element(grad.begin(), grad.end());
  for (uint i = 0; i < grad.size(); i++){
    grad[i] /= max;
  }
}


void construct_linear_gradient(int num_points, vec &grad)
{
  // construct a normalized, centered linear gradient

  // first get the spacing right
  for (int i = 0; i < num_points; i++){
    grad.push_back((double)i);
  }

  // subtract off the average
  double avg = std::accumulate(grad.begin(), grad.end(), 0.0) / grad.size();
  for (uint i = 0; i < grad.size(); i++){
    grad[i] -= avg;
  }

  // normalize by largest value
  double max = *std::max_element(grad.begin(), grad.end());
  for (uint i = 0; i < grad.size(); i++){
    grad[i] /= max;
  }
}


void draw_graph(TGraph gr, string fname, string title)
{
  // Set up the graph
  string new_title(title);
  new_title.append("; time [ms]; amplitude [a.u.]");
  gr.SetTitle(new_title.c_str());

  // Get our own TCanvas
  TCanvas c1;

  // Draw the waveform
  gr.Draw();
  c1.Print(fname.c_str());
}

void draw_graph(const vec &wf, const vec &tm, string fname, string title)
{
  // Create a graph
  TGraph gr(wf.size(), &tm[0], &wf[0]);

  draw_graph(gr, fname, title);
}


void draw_fid(const FID &my_fid, string fname, string title)
{
  // Get the data vectors
  vec wf = my_fid.wf();
  vec tm = my_fid.tm();

  draw_graph(wf, tm, fname, title);
}


void draw_fid_time_fit(const FID &my_fid, string fname, string title)
{
  // Copy the graph
  TGraph gr = my_fid.gr_time_series();
  draw_graph(gr, fname, title);
}


void draw_fid_freq_fit(const FID &my_fid, string fname, string title)
{
  // Copy the graph
  TGraph gr = my_fid.gr_freq_series();
  draw_graph(gr, fname, title);
}

void draw_fid_time_res(const FID &my_fid, string fname, string title)
{
  // Copy the residuals
  vec res = my_fid.res();

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


void draw_fid_freq_res(const FID &my_fid, string fname, string title)
{
  // Copy the residuals
  vec res = my_fid.res();

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

void calc_freq_write_csv(FID& my_fid, ofstream& out)
{
  // Test all the frequency extraction methods and write the results
  out << my_fid.CalcZeroCountFreq() << ", " << 0.0 << ", ";
  out << my_fid.CalcCentroidFreq() << ", " << 0.0 << ", ";
  out << my_fid.CalcAnalyticalFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcLorentzianFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSoftLorentzianFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcExponentialFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSinusoidFreq() << ", " << my_fid.chi2() << endl;;
}

void calc_phase_freq_write_csv(FID& my_fid, ofstream& out)
{
  // Test all the frequency extraction methods using the phase
  out << my_fid.CalcPhaseFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseFreq(2) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseFreq(3) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseDerivFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseDerivFreq(2) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseDerivFreq(3) << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSinusoidFreq() << ", " << my_fid.chi2() << endl;
}

void read_fid_file(string fname, vec &wf, vec &tm)
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

} // ::fid
