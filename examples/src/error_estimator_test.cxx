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

//--- std includes ------------ ----------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
using std::cout;
using std::endl;
using std::string;

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  // Some necessary parameters
  int rounds = 10;

  std::vector<double> wf;
  std::vector<double> tm = time_vector();

  // Ouput filestream.
  std::ofstream out;
  out.precision(12);

  // now get optional output file
  if (argc > 1) {

    out.open(argv[1]);

  } else {

    out.open("data/error_estimator_data.csv");
  }

  // Create random number engine/distribution.
  std::default_random_engine gen;
  std::uniform_real_distribution<double> rand_flat_dist(35.0, 40.0);

  FidFactory ff;
  ff.SetSignalToNoise(100 * 100);

  for (int i = 0; i < rounds; ++i) {
    
    // Make ideal Fid waveform
    ff.SetFidFreq(rand_flat_dist(gen));
    ff.IdealFid(wf, tm);

    Fid myfid(wf, tm);
    out << ff.freq() << "," << 0.0 << ",";
    myfid.WriteMethodCsv(out);
  }

  out.close();

  return 0;
}
