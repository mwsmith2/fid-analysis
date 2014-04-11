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

  fid::ConstructTimeVector(len_fids, i_time, d_time, tm);

  ofstream out;
  out.precision(10);
  out.open("data/ideal_fid_sweep_data.csv");

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

          fid::ideal_fid(wf, tm, f, p, s);
          fid::FID my_fid(wf, tm);

          out << my_fid.CalcZeroCountFreq() << ", ";
          out << my_fid.CalcCentroidFreq() << ", ";
          out << my_fid.CalcAnalyticalFreq() << ", ";
          out << my_fid.chi2() << ", ";
          out << my_fid.CalcLorentzianFreq() << ", ";
          out << my_fid.chi2() << ", ";
          out << my_fid.CalcSoftLorentzianFreq() << ", ";
          out << my_fid.chi2() << ", ";
          out << my_fid.CalcExponentialFreq() << ", ";
          out << my_fid.chi2() << ", ";
          out << my_fid.CalcPhaseFreq() << ", ";
          out << my_fid.chi2() << ", ";
          out << my_fid.CalcSinusoidFreq() << endl;
          out << my_fid.chi2() << ", ";

        } // n_fids

      } // snr

    } // phi

  } // freq

  out.close();
  return 0;
}
