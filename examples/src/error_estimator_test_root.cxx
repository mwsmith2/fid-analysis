/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
file:   error_estimator_test_root.cxx

notes: Tests the techniques used to estimate the fit errors on the
      frequency extraction.

usage:

./error_estimator_test [<output_file>]

The parameters in brackets are optional.  The default output is 
error_estimator_data.csv.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <random>
#include <string>
using std::cout;
using std::endl;
using std::string;
using std::to_string;

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>
#include "TFile.h"
#include "TTree.h"

//--- project includes ------------------------------------------------------//
#include "fid.h"

#define FID_LEN 10000
#define FREQ_METHOD 6

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

  struct fid_data {
    Double_t i_wf;
    Double_t f_wf;
    Double_t i_psd;
    Double_t f_psd;
    Double_t freq_def;
    Double_t wf[FID_LEN];
    Double_t psd[FID_LEN/2];
    Double_t phase[FID_LEN];
    Double_t freq_ext[FREQ_METHOD];
    Double_t freq_err[FREQ_METHOD];
    Double_t fit[FREQ_METHOD][FID_LEN];
  };

  fid_data myfid_data;

  string br_vars;
  br_vars += "i_wf/D:f_wf/D:i_psd/D:f_psd/D:freq_def/D:";
  br_vars += "wf[" + to_string(FID_LEN) + "]/D:";
  br_vars += "psd[" + to_string(FID_LEN/2) + "]/D:";
  br_vars += "phase[" + to_string(FID_LEN) + "]/D:";
  br_vars += "freq_ext[" + to_string(FREQ_METHOD) + "]/D:";
  br_vars += "freq_err[" + to_string(FREQ_METHOD) + "]/D:";
  br_vars += "fit[" + to_string(FREQ_METHOD) + "][" + to_string(FID_LEN) + "]/D";

  TFile pf(out_file.c_str(), "RECREATE");
  TTree pt("t", "fid_fits");
  pt.Branch("fid_data", &myfid_data, br_vars.c_str());

  // some necessary parameters
  std::vector<double> wf;
  std::vector<double> tm;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  tm = construct_range(sim::start_time, final_time, sim::delta_time);

  // Create random number engine/distribution.
  std::default_random_engine gen;
  std::uniform_real_distribution<double> rand_flat_dist(35.0, 40.0);

  for (int i = 0; i < 1000; ++i) {
    
    if (i % 250 == 0) {
      cout << "Processing round " << i << "." << endl;
    }

    // Make ideal FID waveform
    double freq = rand_flat_dist(gen);
    ideal_fid(wf, tm, freq, 0.0, 1000);

    FID myfid(wf, tm);
    myfid_data.i_wf = myfid.i_wf();
    myfid_data.f_wf = myfid.f_wf();
    myfid_data.i_psd = myfid.i_fft();
    myfid_data.f_psd = myfid.f_fft();
    myfid_data.freq_def = freq;
    std::copy(myfid.wf().begin(), myfid.wf().end(), myfid_data.wf);
    std::copy(myfid.power().begin(), myfid.power().end(), myfid_data.psd);
    std::copy(myfid.phase().begin(), myfid.phase().end(), myfid_data.phase);

    int idx = 0;
    myfid_data.freq_ext[idx] = myfid.CalcZeroCountFreq();
    myfid_data.freq_err[idx] = myfid.freq_err();

    idx++;
    myfid_data.freq_ext[idx] = myfid.CalcCentroidFreq();
    myfid_data.freq_err[idx] = myfid.freq_err();

    idx++;
    myfid_data.freq_ext[idx] = myfid.CalcAnalyticalFreq();
    myfid_data.freq_err[idx] = myfid.freq_err();

    for (int i = myfid.i_fft(); i < myfid.f_fft(); ++i) {
      myfid_data.fit[idx][i] = myfid.f_fit().Eval(myfid.fftfreq()[i]);
    }

    idx++;
    myfid_data.freq_ext[idx] = myfid.CalcLorentzianFreq();
    myfid_data.freq_err[idx] = myfid.freq_err();

    for (int i = myfid.i_fft(); i < myfid.f_fft(); ++i) {
      myfid_data.fit[idx][i] = myfid.f_fit().Eval(myfid.fftfreq()[i]);
    }

    idx++;
    myfid_data.freq_ext[idx] = myfid.CalcPhaseFreq();
    myfid_data.freq_err[idx] = myfid.freq_err();

    for (int i = myfid.i_wf(); i < myfid.f_wf(); ++i) {
      myfid_data.fit[idx][i] = myfid.f_fit().Eval(myfid.tm()[i]);
    }

    idx++;
    myfid_data.freq_ext[idx] = myfid.CalcSinusoidFreq();
    myfid_data.freq_err[idx] = myfid.freq_err();

    for (int i = myfid.i_wf(); i < myfid.f_wf(); ++i) {
      myfid_data.fit[idx][i] = myfid.f_fit().Eval(myfid.tm()[i]);
    }

    pt.Fill();
  }

  pf.Write();
  pf.Close();  
  return 0;
}
