/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
date:   2015/01/09
file:   single_fid_analyzer.cxx

notes: This program extracts the frequency of a single fid file and 
      produces plots of the Fid, envelope function and fits. The analysis
      is extensive testing all available methods.

usage:

./single_fid_analyzer <fid_data> [<output_dir>]

The parameters in brackets are optional.  The default output is 
data/ and the directory will be created if not present.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <cassert>
using std::string;

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  // Make sure a data file was specified and get the file handle.
  if (argc < 2) {
    std::cout << "Insufficient arguments: input data must be specified.";
    std::cout << std::endl;
    std::cout << "Usage: ./single_fid_analyzer <input-file>" << std::endl;
    exit(1);
  }

  // Load the data.
  string data_file(argv[1]);

  // Now check for an optional output directory.
  string out_dir;

  if (argc > 2){

    out_dir = string(argv[2]);

  } else {

    out_dir = string("data/");
  }

  // Create an output file handle.
  string out_file(out_dir + string("single_fid_data.csv"));

  // Make a string to hold the location of the figure output directory.
  string fig_dir(out_dir + string("fig/"));

  // Make certain that the directory exists by creating it.
  boost::filesystem::path dir(fig_dir);
  boost::filesystem::create_directories(dir);

  // Read the Fid data and create Fid object.
  Fid my_fid(data_file);

  // Open the output filestream.
  std::ofstream out;
  out.precision(10);
  out.open(out_file);

  // Test all methods and save the frequency results
  my_fid.WriteFreqCsv(out);

  // Make plots of the fits and residuals.
  string title("Fid Fit");
  my_fid.GetFreq("analytical");
  my_fid.SaveFreqFit(fig_dir + string("analytical_fit.pdf"), title);
  my_fid.SaveFreqRes(fig_dir + string("analytical_res.pdf"), title);  

  my_fid.GetFreq("lorentzian");
  my_fid.SaveFreqFit(fig_dir + string("lorentzian_fit.pdf"), title);
  my_fid.SaveFreqRes(fig_dir + string("lorentzian_res.pdf"), title);  

  my_fid.GetFreq("exponential");
  my_fid.SaveFreqFit(fig_dir + string("exponential_fit.pdf"), title);
  my_fid.SaveFreqRes(fig_dir + string("exponential_res.pdf"), title);  

  my_fid.GetFreq("phase");
  my_fid.SaveTimeFit(fig_dir + string("lin_phase_fit.pdf"), title);
  my_fid.SaveTimeRes(fig_dir + string("lin_phase_res.pdf"), title);  

  my_fid.GetFreq("sinusoid");
  my_fid.SaveTimeFit(fig_dir + string("sinusoid_fit.pdf"), title);
  my_fid.SaveTimeRes(fig_dir + string("sinusoid_res.pdf"), title);  

  // Close the output filestream
  out.close();
  return 0;
}
