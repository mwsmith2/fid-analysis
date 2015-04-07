#ifndef FID_INCLUDE_FID_UTILS_H_
#define FID_INCLUDE_FID_UTILS_H_

/*===========================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes: 

	This library consists of several frequency extraction and analysis 
	methods for FIDs as well as a class to encapsulate all the ideas.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
using std::vector;
using std::cout;
using std::endl;

//--- other includes --------------------------------------------------------//
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"

//--- project includes ------------------------------------------------------//
#include "fid_class.h"
#include "fid_sim.h"
#include "fid_math.h"

namespace fid
{
	// Classes defined in separate headers
	class FidFactory;
	class FID;

  // Run this function first thing in any module to load a custom configuration
  void load_params(string conf_file);

	// Declare utility functions.
	void ideal_fid(vec& wf, vec& tm, double f, double phi=0.0, double snr=100.0, double tau=10.0, double t0=0.0);

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
  void calc_freq_write_csv(FID &my_fid, ofstream &out);

  // Try all phase related frequency extraction methods and write to csv
  void calc_phase_freq_write_csv(FID &my_fid, ofstream &out);

  // Read a FID file which is two space delimited columns (time voltage)
  void read_fid_file(string fname, vec &wf, vec &tm);

  // Get a time vector for the FID
  void construct_time_vector(int num_times, double t0, double dt, vec &tm);

  // Get a linear gradient with a max(abs(grad)) == 1 and mean(grad) == 0
  void construct_linear_gradient(int num_points, vec &grad);

  // Get a quadratic gradient with a max(abs(grad)) == 1 and mean(grad) = 0
  void construct_quadratic_gradient(int num_points, vec &grad);
  
} // fid

#endif