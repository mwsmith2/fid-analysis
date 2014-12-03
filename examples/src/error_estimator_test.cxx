/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
date:   2014/12/02

notes: Tests the techniques used to estimate the fit errors on the
      frequency extraction.

usage:

./error_estimator_test [<output_file>]

The parameters in brackets are optional.  The default output is 
error_estimator_data.csv.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <random>

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  // filenames
  string out_file;
  string out_dir;
  string fig_dir;

  // now get optional output file
  if (argc > 1) {

    out_file = string(argv[1]);

  } else {

    out_file = string("error_estimator_data.root");

  }

  // open the output file
  ofstream out;
  out.precision(12);
  out.open(out_file);

  // some necessary parameters
  vec wf;
  vec tm;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  tm = construct_range(sim::start_time, final_time, sim::delta_time);

  // Create random number engine/distribution.
  std::default_random_engine gen;
  std::uniform_real_distribution<double> rand_flat_dist(35.0, 40.0);

  for (int i = 0; i < 2000; ++i) {
    
    // Make ideal FID waveform
    double freq = rand_flat_dist(gen);
    ideal_fid(wf, tm, freq, 0.0, 100 * 100);

    FID myfid(wf, tm);
    out << freq << "\t" << 0.0 << "\t";
    out << myfid.CalcZeroCountFreq() << "\t" << myfid.freq_err() << "\t";
    out << myfid.CalcCentroidFreq() << "\t" << myfid.freq_err() << "\t";
    out << myfid.CalcAnalyticalFreq() << "\t" << myfid.freq_err() << "\t";
    out << myfid.CalcLorentzianFreq() << "\t" << myfid.freq_err() << "\t";
    out << myfid.CalcPhaseFreq() << "\t" << myfid.freq_err() << "\t";
    out << myfid.CalcSinusoidFreq() << "\t" << myfid.freq_err() << endl;;

}
  // todo residuals
  out.close();
  return 0;
}
