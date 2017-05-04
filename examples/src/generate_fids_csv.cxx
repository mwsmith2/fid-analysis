/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
file:   generate_fids_csv.cxx

about: This is a new test program for my Fid libraries 

\*===========================================================================*/


//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//--- other includes --------------------------------------------------------//

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;
using std::cout;
using std::endl;


int main(int argc, char **argv)
{
  // initialize the configurable parameters
  if (argc > 1) load_params(argv[1]);

  int fid_idx = 0;
  char outfile[256];

  // some necessary parameters
  std::vector<double> wf;
  std::vector<double> tm = time_vector();
  std::vector<double> freqs;
  std::vector<double> snrs;

  // Set the ranges.
  freqs = construct_range(23.0, 24.0, 0.1);
  snrs = construct_range(10000.0, 100000.0, 30000.0);

  // Make FidFactory
  FidFactory ff;

  // begin sweeps
  for (auto f: freqs) {

    for (auto s: snrs) {

      if (freqs.size() > 1) cout << "Running for frequency " << f;
      if (snrs.size() > 1) cout << ", signal-to-noise " << s;
      cout << endl;

      ff.SimulateFid(wf, tm);

      Fid my_fid(wf, tm);

      snprintf(outfile, 256, "data/sim_%05i.fid", fid_idx++);
      my_fid.SaveData(outfile);

    } // snr

  } // freq

  return 0;
}
