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

  // Replicate the functionality of MATLAB/numpy linspace
  std::vector<double> linspace(double x0, double xf, int n=0);

  // Get a linear gradient with a max(abs(grad)) == 1 and mean(grad) == 0
  std::vector<double> normalized_gradient(int npoints, int poln=1); 

} // fid

#endif