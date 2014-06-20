/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
date:   15/04/14

notes: This program extracts the frequency of a single fid file and produces plots of the FID, envelope function and fits.

usage:

./single_fid_analyzer <fid_data> [<output_dir>]

The parameters in brackets are optional.  The default output is 
data/single_fid.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <cassert>

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>

//--- project includes ------------------------------------------------------//
#include "fid.h"
#include "fid_class.h"
#include "fid_params.h"
#include "fid_utilities.h"

using namespace fid;
using namespace fid::sweep;

int main(int argc, char **argv)
{
  // initialize the configurable parameters
  load_params(argc, argv);

  // filenames
  string data_file;
  string out_file;
  string out_dir;
  string fig_dir;

  // get fid data file
  assert(argc > 1);
  data_file = string(argv[1]);

  // now get optional output file
  if (argc > 2){

    out_dir = string(argv[2]);

  } else {

    out_dir = string("data/single_fid/");

  }

  out_file = out_dir + string("single_fid_data.csv");
  fig_dir = out_dir + string("fig/");

  // open the output file  
  ofstream out;
  out.precision(10);
  out.open(out_file);

  // make certain that the directory exists by creating it
  boost::filesystem::path dir(fig_dir);
  boost::filesystem::create_directories(dir);

  // some necessary parameters
  vec wf;
  vec tm;

  // read the FID data and create FID object
  read_fid_file(data_file, wf, tm);
  FID my_fid(wf, tm);

  // test all methods and save the frequency results
  calc_freq_write_csv(my_fid, out);

  out.close();

  // now make some fit plots
  string title("FID Fit");
  my_fid.CalcAnalyticalFreq();
  draw_fid_freq_fit(my_fid, fig_dir + string("analytical_fit.pdf"), title);
  draw_fid_freq_res(my_fid, fig_dir + string("analytical_res.pdf"), title);  

  my_fid.CalcLorentzianFreq();
  draw_fid_freq_fit(my_fid, fig_dir + string("lorentzian_fit.pdf"), title);
  draw_fid_freq_res(my_fid, fig_dir + string("lorentzian_res.pdf"), title);  

  my_fid.CalcExponentialFreq();
  draw_fid_freq_fit(my_fid, fig_dir + string("exponential_fit.pdf"), title);
  draw_fid_freq_res(my_fid, fig_dir + string("exponential_res.pdf"), title);  

  my_fid.CalcPhaseFreq();
  draw_fid_time_fit(my_fid, fig_dir + string("lin_phase_fit.pdf"), title);
  draw_fid_time_res(my_fid, fig_dir + string("lin_phase_res.pdf"), title);  

  my_fid.CalcSinusoidFreq();
  draw_fid_time_fit(my_fid, fig_dir + string("sinusoid_fit.pdf"), title);
  draw_fid_time_res(my_fid, fig_dir + string("sinusoid_res.pdf"), title);  

  // todo residuals
  out.close();
  return 0;
}