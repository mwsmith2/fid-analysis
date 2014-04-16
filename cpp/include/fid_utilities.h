#ifndef FID_ANALYSIS_INCLUDE_FID_UTILITIES_H_
#define FID_ANALYSIS_INCLUDE_FID_UTILITIES_H_



//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
using std::vector;
using std::string;
using std::ofstream;

//--- root includes ---------------------------------------------------------//
#include "TCanvas.h"
#include "TGraph.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_class.h"
#include "fid.h"

namespace fid{

  // Plots and saves an image
  void draw_graph(TGraph gr, string fname, string title);
  void draw_graph(const vec &wf, const vec &tm, string fname, string title);

  // Plots and save specific FID images
  void draw_fid(const FID &my_fid, string fname, string title);
  void draw_fid_time_fit(const FID &my_fid, string fname, string title);
  void draw_fid_freq_fit(const FID &my_fid, string fname, string title);
  void draw_fid_time_res(const FID &my_fid, string fname, string title);
  void draw_fid_freq_res(const FID &my_fid, string fname, string title);

  // Try all frequency extraction methods and write in a csv format
  void calc_freq_save_csv(FID &my_fid, ofstream &out);

  // Read a FID file which is two space delimited columns (time voltage)
  void read_fid_file(string fname, vec &wf, vec &tm);

  // Add Gaussian noise to the given waveform
  void add_white_noise(vec &wf, double snr=100.0);

  // Get a time vector for the FID
  void construct_time_vector(int num_times, double t0, double dt, vec &tm);

  // Get a linear gradient with a max(abs(grad)) == 1 and mean(grad) == 0
  void construct_linear_gradient(int num_points, vec &grad);

  // Get a quadratic gradient with a max(abs(grad)) == 1 and mean(grad) = 0
  void construct_quadratic_gradient(int num_points, vec &grad);

  // Use for getting sweep ranges
  template<typename T>
  inline vector<T> construct_sweep_range(const vector<T> &range_vec){
    vector<T> res;

    for (T it = range_vec[0]; it <= range_vec[1]; it += range_vec[2]){
      res.push_back(it);
    }

    return res;
  }

} // fid

#endif