/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
date:   15/04/14

notes: This program extracts the frequency of a single fid file and produces plots of the FID, envelope function and fits.

usage:

./single_fid_analyzer <fid_data> [<output_file> <figure_directory>]

The parameters in brackets are optional.

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
  string fig_dir;

  // get fid data file
  assert(argc > 1);
  data_file = string(argv[1]);

  // now get optional output file
  if (argc > 2){

    out_file = string(argv[2]);

  } else {

    out_file = string("data/ideal_fid_sweep_data.csv");

  }

  ofstream out;
  out.precision(10);
  out.open(out_file);

  // now get optional figure directory
  if (argc > 3){

    fig_dir = string(argv[3]);

  } else {

    fig_dir = string("data/fig/single_fid/");

  }

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
  calc_freq_save_csv(my_fid, out);

  // now make some plots
  string title("FID Fit");
  my_fid.CalcAnalyticalFreq();
  draw_fid_freq_fit(my_fid, fig_dir + string("/analytical_fit.pdf"), title);
  draw_fid_freq_fit(my_fid, fig_dir + string("/lorentzian_fit.pdf"), title);
  draw_fid_freq_fit(my_fid, fig_dir + string("/exponential_fit.pdf"), title);
  draw_fid_time_fit(my_fid, fig_dir + string("/lin_phase_fit.pdf"), title);
  draw_fid_time_fit(my_fid, fig_dir + string("/sin_fit.pdf"), title);

  // todo residuals
  out.close();
  return 0;
}
