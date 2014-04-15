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

  // Plots and saves an image of the FID
  void draw_fid(const vec &wf, const vec &tm, 
    const string filename, const string title);
  void draw_fid(FID &my_fid, const string filename, const string title);

  void calc_freq_save_csv(FID &my_fid, ofstream &out);

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