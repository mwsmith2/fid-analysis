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

//--- other includes --------------------------------------------------------//
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

//--- project includes ------------------------------------------------------//
#include "fid.h"

using std::cout;
using std::endl;
using namespace fid;


int main(int argc, char **argv)
{
  // Initialize the configurable parameters.
  if (argc > 1) load_params(argv[1]);

  // Standard variables.
  std::vector<double> wf(0.0, sim::num_samples);
  std::vector<double> tm(0.0, sim::num_samples);

  int num_fids = 10;
  double max_grad = 500; // in ppb
  double delta_b = 0.0;
  double center_b = 47.0;

  double final_time = sim::start_time + sim::num_samples * sim::sample_time;
  tm = construct_range(sim::start_time, final_time, sim::sample_time);

  // Set up the ROOT tree to hold the results
  TFile pf("sim_fids.root", "recreate");
  TTree pt("t", "Fid Tree");
  cout.precision(12);

  pt.Branch("db", &delta_b, "delta_b/D");
  pt.Branch("fid", &wf[0], TString::Format("fid[%d]/D", sim::num_samples));
  pt.Branch("tm", &tm[0], TString::Format("time[%d]/D", sim::num_samples));

  FidFactory ff;

  // Begin simulating FIDs.
  cout << "num_fids: " << num_fids << endl;

  for (int i = -1 * num_fids / 2; i < num_fids / 2 + 1; ++i) {

  	// Set the FID frequency.
    delta_b = max_grad * 2.0 * i / num_fids;

    ff.SetFidFreq((1.0 + 1.0e-9 * delta_b) * center_b);
    ff.SimulateFid(wf, tm);

    pt.Fill();

  } // grad

  pf.Write();
  pf.Close();

  return 0;
}
