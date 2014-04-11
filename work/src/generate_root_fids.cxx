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

  fid::ConstructTimeVector(len_fids, i_time, d_time, tm);

  vec freqs;
  vec phases;
  vec snrs;

  // Get the range to sweep over.
  if (freq_sweep){
    freqs = ConstructSweepRange(freq_range);
  } else {
    freqs.push_back(s_freq);
  }

  if (phase_sweep){
    phases = ConstructSweepRange(phase_range);
  } else {
    phases.push_back(s_phase);
  }

  if (snr_sweep){
    snrs = ConstructSweepRange(snr_range);
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

        cout << "SimulatedFid: ";
        for (auto it = wf.begin(); it != wf.end(); ++it){
          cout << *it;
        }
        cout << endl;
        fid::FID my_fid(wf, tm);

      } // snr

    } // phi

  } // freq

  return 0;
}
