#include "fid_utilities.h"

namespace fid {

void construct_time_vector(int num_times, double t0, double dt, vec &tm)
{
  if (tm.size() != num_times){
    tm.resize(num_times);
  }

  for (int i = 0; i < num_times; i++){
    tm[i] = dt * i + t0;
  }
}


void construct_quadratic_gradient(int num_points, vec &grad)
{
  // construct a normalize linear gradient

  // first get the spacing right
  for (int i = 0; i < num_points; i++){
    grad.push_back((double)i * i);
  }

  // subtract off the average
  double avg = std::accumulate(grad.begin(), grad.end(), 0.0) / grad.size();
  for (int i = 0; i < grad.size(); i++){
    grad[i] -= avg;
  }

  // normalize by largest value
  double max = *std::max_element(grad.begin(), grad.end());
  for (int i = 0; i < grad.size(); i++){
    grad[i] /= max;
  }
}


void construct_linear_gradient(int num_points, vec &grad)
{
  // construct a normalize linear gradient

  // first get the spacing right
  for (int i = 0; i < num_points; i++){
    grad.push_back((double)i * i);
  }

  // subtract off the average
  double avg = std::accumulate(grad.begin(), grad.end(), 0.0) / grad.size();
  for (int i = 0; i < grad.size(); i++){
    grad[i] -= avg;
  }

  // normalize by largest value
  double max = *std::max_element(grad.begin(), grad.end());
  for (int i = 0; i < grad.size(); i++){
    grad[i] /= max;
  }
}


void draw_fid(const vec &wf, 
                  const vec &tm, 
                  const string filename,
                  const string title)
{
  // Get our own TCanvas
  TCanvas c1;

  // Set up the graph
  TGraph gr(wf.size(), &tm[0], &wf[0]);
  string new_title(title);
  new_title.append("; time [ms]; amplitude [a.u.]");
  gr.SetTitle(title.c_str());

  // Draw the waveform
  gr.Draw();
  c1.Print(filename.c_str());
}


void draw_fid(FID &my_fid, const string filename, const string title)
{
  // Get our own TCanvas
  TCanvas c1;

  // Set up the graph
  int N = my_fid.wf().size();
  TGraph gr(N, &my_fid.tm()[0], &my_fid.wf()[0]);
  string new_title(title);
  new_title.append("; time [ms]; amplitude [a.u.]");
  gr.SetTitle(title.c_str());

  // Draw the waveform
  gr.Draw();
  c1.Print(filename.c_str());
}


void add_white_noise(vec &wf, double snr){
  static std::default_random_engine gen;
  static std::normal_distribution<double> nrm(0.0, snr);

  double max = *std::max_element(wf.begin(), wf.end());
  double min = *std::min_element(wf.begin(), wf.end());
  double scale = max > min ? max : min;

  for (auto x : wf){
    x += scale * nrm(gen);
  }
}

void calc_freq_save_csv(FID& my_fid, ofstream& out)
{
  // Test all the frequency extraction methods and write the results
  out << my_fid.CalcZeroCountFreq() << ", " << 0.0 << ", ";
  out << my_fid.CalcCentroidFreq() << ", " << 0.0 << ", ";
  out << my_fid.CalcAnalyticalFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcLorentzianFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSoftLorentzianFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcExponentialFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcPhaseFreq() << ", " << my_fid.chi2() << ", ";
  out << my_fid.CalcSinusoidFreq() << ", " << my_fid.chi2() << ", ";
}

} // fid





