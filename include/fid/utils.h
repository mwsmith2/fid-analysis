#ifndef FID_INCLUDE_FID_UTILS_H_
#define FID_INCLUDE_FID_UTILS_H_

/*==========================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes: 

	This library consists of several frequency extraction and analysis 
	methods for FIDs as well as a class to encapsulate all the ideas.

\*==========================================================================*/

//--- std includes ---------------------------------------------------------//
#include <string>
#include <vector>
#include <random>

//--- other includes -------------------------------------------------------//
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "TGraph.h"
#include "TCanvas.h"

//--- project includes -----------------------------------------------------//
#include "fid/fid.h"

namespace fid
{
	// Classes ared defined in a separate header.
	class FidFactory;
	class FID;

  // Declare utility functions.

  // Run this function first thing to load a custom configuration.
  void load_params(std::string conf_file);

  // Plots and saves the image.
  void draw_graph(TGraph gr, std::string fname, std::string title);

  // Plots and saves an image of the FID.
  void draw_graph(const std::vector<double> &wf, 
                  const std::vector<double> &tm, 
                  std::string fname, 
                  std::string title);

  // Plots and save specific FID images
  void draw_fid(const FID &my_fid, std::string fname, std::string title);

  void draw_fid_time_fit(const FID &my_fid, 
                         std::string fname, 
                         std::string title);

  void draw_fid_freq_fit(const FID &my_fid, 
                         std::string fname, 
                         std::string title);

  void draw_fid_time_res(const FID &my_fid, 
                         std::string fname, 
                         std::string title);

  void draw_fid_freq_res(const FID &my_fid, 
                         std::string fname, 
                         std::string title);


  // Try all frequency extraction methods and write in a csv format
  void calc_freq_write_csv(FID &my_fid, std::ofstream &out);

  // Try all phase related frequency extraction methods and write to csv
  void calc_phase_freq_write_csv(FID &my_fid, std::ofstream &out);

  // Read a FID file which is two space delimited columns (time voltage)
  void read_fid_file(std::string fname, 
                     std::vector<double> &wf, 
                     std::vector<double> &tm);

  // Read a FID file which is two space delimited columns (time voltage)
  void write_fid_file(std::string fname, 
                      const std::vector<double> &wf, 
                      const std::vector<double> &tm);


  // Get a time vector for the FID
  void construct_time_vector(int num_times, 
                             double t0, 
                             double dt, 
                             std::vector<double> &tm);

  // Get a linear gradient with a max(abs(grad)) == 1 and mean(grad) == 0
  void construct_linear_gradient(int noints, std::vector<double> &grad);

  // Get a quadratic gradient with a max(abs(grad)) == 1 and mean(grad) = 0
  void construct_quadratic_gradient(int npoints, std::vector<double> &grad);
  
} // fid

#endif