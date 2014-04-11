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

namespace fid {

int main(int argc, char **argv)
{
  using namespace fid::sweep;

  // initialize the configurable parameters
  load_params(argc, argv);

  // some necessary parameters
  vec wf;
  vec tm;

  ofstream out;
  out.precision(10);

  fid::ConstructTimeVector(len_fids, i_time, d_time, tm);

  out.open("ideal_fid_sweep_data.csv");

  bool freq_sweep = true;
  bool phi_sweep = false;
  bool snr_sweep = false;

  // fix loop starts
  i_freq = freq_sweep ? i_freq : freq;
  i_phi  = phi_sweep ? i_phi : phi;
  i_snr  = snr_sweep ? i_snr : snr;      

  // fix loop maxes
  f_freq = freq_sweep ? f_freq : freq;
  f_phi  = phi_sweep ? f_phi : phi;
  f_snr  = snr_sweep ? f_snr : snr;

  // begin sweeps
  for (double freq = i_freq; freq <= f_freq; freq += d_freq){

    if (freq_sweep) cout << "Running for frequency " << freq;
    out << freq << ", ";

    for (double phi = i_phi; phi <= f_phi; phi += d_phi){

      if (phi_sweep) cout << ", phase " << phi;
      out << phi << ", ";

      for (double snr = i_snr; snr <= f_snr; snr += d_snr){

        if (snr_sweep) cout << ", signal-to-noise " << snr;
        out << snr << ", ";

        cout  << ".\n";

        for (int i = 0; i < num_fids; i++){

            fid::ideal_fid(wf, tm, freq, phi, snr);

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

} // ::fid