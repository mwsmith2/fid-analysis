/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This is a new test program for my FID libraries 

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>

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

  // some necessary parameters
  vec wf;
  vec tm;

  fid::construct_time_vector(len_fids, i_time, d_time, tm);

  ofstream out;
  out.precision(10);
  out.open("data/ideal_fid_sweep_data.csv");

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

  // begin sweeps
  for (auto f: freqs){

    for (auto p: phases){

      for (auto s: snrs){

        if (freq_sweep) cout << "Running for frequency " << f;
        if (phase_sweep) cout << ", phase " << p;
        if (snr_sweep) cout << ", signal-to-noise " << s;
        cout << endl;

        for (int i = 0; i < num_fids; i++){

          if (freq_sweep) out << f << ", ";
          if (phase_sweep) out << p << ", ";
          if (snr_sweep) out << s << ", ";

          ideal_fid(wf, tm, f, p, s);
          FID my_fid(wf, tm);

          calc_freq_save_csv(my_fid, out);
        } // n_fids

      } // snr

    } // phi

  } // freq

  out.close();
  return 0;
}
