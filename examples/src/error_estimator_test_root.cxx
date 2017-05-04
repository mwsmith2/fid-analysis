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

//--- other includes --------------------------------------------------------//
#include <boost/filesystem.hpp>
#include "TFile.h"
#include "TTree.h"

//--- project includes ------------------------------------------------------//
#include "fid.h"

#define FID_LEN 10000
#define FREQ_METHOD 6

using std::cout;
using std::endl;
using std::string;
using std::to_string;
using namespace fid;


int main(int argc, char **argv)
{
  // Necessary variables.
  int rounds = 10;
  std::vector<double> wf;
  std::vector<double> tm = time_vector(); // default time vector in fid::

  // Ouput filename.
  string outfile;

  // now get optional output file
  if (argc > 1) {

    outfile = string(argv[1]);

  } else {

    outfile = string("data/error_estimator_data.root");

  }

  fid_freq_t fid_data;
  Double_t freq;

  TFile pf(outfile.c_str(), "RECREATE");
  TTree pt("t", "FID Fits");
  pt.Branch("data", &fid_data, fid_freq_str);
  pt.Branch("seed", &freq, "freq/D");

  // Create random number engine/distribution.
  std::default_random_engine gen;
  std::uniform_real_distribution<double> rand_flat_dist(35.0, 40.0);

  FidFactory ff;
  ff.SetSignalToNoise(100*100);
  ff.SetMixdownPhi(0.0);

  for (int i = 0; i < rounds; ++i) {
    
    if (i % (rounds / 5) == 0) {
      cout << "Processing round " << i << " of " << rounds << "." << endl;
    }

    // Make ideal Fid waveform
    freq = rand_flat_dist(gen);

    ff.SetFidFreq(freq);
    ff.IdealFid(wf, tm);
    Fid myfid(wf, tm);

    myfid.CopyStruct(fid_data);

    pt.Fill();
  }

  pf.Write();
  pf.Close();  
  return 0;
}
