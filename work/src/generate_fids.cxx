/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This is a new test program for my FID libraries 

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>

//--- other includes --------------------------------------------------------//
#include "TFile.h"
#include "TTree.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_class.h"
#include "fid_utilities.h"
#include "fid_sim.h"

using namespace fid;
using namespace fid::sweep;

int main(int argc, char **argv)
{
  // initialize the configurable parameters
  load_params(argc, argv);

  // some necessary parameters
  vec wf;
  vec tm;

  construct_time_vector(len_fids, i_time, d_time, tm);

  vec freqs;
  vec phases;
  vec snrs;

  // Get the range to sweep over.
  if (freq_sweep){
    freqs = construct_sweep_range(freq_range);
  } else {
    freqs.push_back(s_freq);
  }

  if (phase_sweep){
    phases = construct_sweep_range(phase_range);
  } else {
    phases.push_back(s_phase);
  }

  if (snr_sweep){
    snrs = construct_sweep_range(snr_range);
  } else {
    snrs.push_back(s_snr);
  }

  // Make FidFactory
  FidFactory ff;

  // Set up the ROOT tree to hold the results
  // todo

  // begin sweeps
  for (auto f: freqs){

    for (auto p: phases){

      for (auto s: snrs){

        if (freq_sweep) cout << "Running for frequency " << f;
        if (phase_sweep) cout << ", phase " << p;
        if (snr_sweep) cout << ", signal-to-noise " << s;
        cout << endl;

        ff.SimulateFid(wf, tm);

        FID my_fid(wf, tm);
        draw_fid(my_fid, string("data/fig/test.pdf"), string("Test FID"));

      } // snr

    } // phi

  } // freq

  return 0;
}
