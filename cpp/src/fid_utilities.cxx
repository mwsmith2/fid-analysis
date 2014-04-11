#include "fid_utilities.h"

namespace fid {

  void ConstructTimeVector(int num_times, double t0, double dt, vec &tm)
  {
    if (tm.size() != num_times){
      tm.resize(num_times);
    }

    for (int i = 0; i < num_times; i++){
      tm[i] = dt * i + t0;
    }

    return;
  }


  void ConstructQuadraticGradient(int num_points, vec &grad)
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

    return;
  }


  void ConstructLinearGradient(int num_points, vec &grad)
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

    return;
  }


  void DrawFID(const vec &wf, 
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

    return;
  }


  void DrawFID(FID &my_fid, const string filename, const string title)
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

    return;
  }


  void AddWhiteNoise(vec &wf, double snr){
    static std::default_random_engine gen;
    static std::normal_distribution<double> nrm(0.0, snr);

    double max = *std::max_element(wf.begin(), wf.end());
    double min = *std::min_element(wf.begin(), wf.end());
    double scale = max > min ? max : min;

    for (auto x : wf){
      x += scale * nrm(gen);
    }

    return;
  }

} // fid





