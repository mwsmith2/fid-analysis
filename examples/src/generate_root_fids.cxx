/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
file:   generate_root_fids.cxx

notes: This is a new test program for my Fid libraries 

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <iostream>
#include <vector>
using std::cout;
using std::endl;

//--- other includes --------------------------------------------------------//
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  // initialize the configurable parameters
  if (argc > 1) load_params(argv[1]);

  // some necessary parameters
  std::vector<double> wf(0.0, sim::num_samples);
  std::vector<double> tm = time_vector();

  double delta_b;
  double max_grad = 500; // in ppb
  int num_fids = 100;
  double freq_0 = sim::larmor_freq;

  // Set up the ROOT tree to hold the results
  TFile pf("sim_fids.root", "recreate");
  TTree pt("t", "Fid Tree");
  cout.precision(12);

  pt.Branch("db", &delta_b, "delta_b/D");
  pt.Branch("fid", &wf[0], TString::Format("fid[%d]/D", sim::num_samples));
  pt.Branch("tm", &tm[0], TString::Format("time[%d]/D", sim::num_samples));

  // begin grad sims
  cout << "num_fids: " << num_fids << endl;
  for (int i = -1 * num_fids / 2; i < num_fids / 2 + 1; ++i){

    delta_b = max_grad * 2.0 * i / num_fids;
    sim::larmor_freq = (1.0 + 1.0e-9 * delta_b) * freq_0;

    std::cout << sim::larmor_freq << std::endl;

    // Make FidFactory
    FidFactory ff;
    ff.SimulateFid(wf, tm);

    pt.Fill();

  } // grad

  pf.Write();
  pf.Close();

  return 0;
}
