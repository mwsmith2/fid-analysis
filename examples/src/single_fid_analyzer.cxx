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
#include <iostream>

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using std::string;
using std::cout;
using std::endl;
using namespace fid;

int main(int argc, char **argv)
{
  string datafile;
  string outdir;
  string outfile;
  string figdir;

  std::ofstream out;
  out.precision(10);

  // Make sure a data file was specified and get the file handle.
  if (argc < 2) {
    cout << "Insufficient arguments: input data must be specified." << endl;
    cout << "Usage: ./single_fid_analyzer <input-file> [output-dir]" << endl;
    exit(1);
  }

  // Load the data.
  datafile = string(argv[1]);

  // Now check for an optional output directory.
  if (argc > 2){

    outdir = string(argv[2]);

  } else {

    outdir = string("data/");
  }

  // Create an output file handle.
  outfile = outdir + string("single_fid_data.csv");

  // Make a string to hold the location of the figure output directory.
  figdir = outdir + string("fig/");

  // Make certain that the directory exists by creating it.
  boost::filesystem::path dir(figdir);
  boost::filesystem::create_directories(dir);

  // Read the Fid data and create Fid object.
  Fid my_fid(datafile);

  // Open the output filestream.
  out.open(outfile);
  my_fid.WriteFreqCsv(out);

  // Make plots of the fits and residuals.
  string title("Fid Fit");
  my_fid.GetFreq("analytical");
  my_fid.SaveFreqFit(dir.string() + string("analytical_fit.png"), title);
  my_fid.SaveFreqRes(dir.string() + string("analytical_res.png"), title);  

  my_fid.GetFreq("lorentzian");
  my_fid.SaveFreqFit(dir.string() + string("lorentzian_fit.png"), title);
  my_fid.SaveFreqRes(dir.string() + string("lorentzian_res.png"), title);  

  my_fid.GetFreq("exponential");
  my_fid.SaveFreqFit(dir.string() + string("exponential_fit.png"), title);
  my_fid.SaveFreqRes(dir.string() + string("exponential_res.png"), title);  

  my_fid.GetFreq("phase");
  my_fid.SaveTimeFit(dir.string() + string("lin_phase_fit.png"), title);
  my_fid.SaveTimeRes(dir.string() + string("lin_phase_res.png"), title);  

  my_fid.GetFreq("sinusoid");
  my_fid.SaveTimeFit(dir.string() + string("sinusoid_fit.png"), title);
  my_fid.SaveTimeRes(dir.string() + string("sinusoid_res.png"), title);  

  // Close the output filestream
  out.close();

  return 0;
}
