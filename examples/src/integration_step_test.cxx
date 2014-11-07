/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
file:   integration_step_test.cxx

about: This is a new test program for my FID libraries 

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <iostream>

//--- other includes --------------------------------------------------------//
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

//--- project includes ------------------------------------------------------//
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  string root_file = "dt_test_fids.root";

  // initialize the configurable parameters
  if (argc > 1) load_params(argv[1]);

  // some necessary parameters
  vec wf;
  vec tm;
  wf.reserve(sim::num_samples);
  tm.reserve(sim::num_samples);
  double dt;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  tm = construct_range(sim::start_time, final_time, sim::delta_time);

  // Set up the ROOT tree to hold the results
  TFile pf(root_file.c_str(), "recreate");
  TTree pt("t", "FID Tree");
  cout.precision(12);

  pt.Branch("dt", &dt, "time_step/D");
  pt.Branch("fid", &wf[0], TString::Format("fid[%d]/D", sim::num_samples));
  pt.Branch("tm", &tm[0], TString::Format("time[%d]/D", sim::num_samples));

  vec stepsizes;
  for (int i = 3; i <= 6; ++i) {
    stepsizes.push_back(2.0 * pow(10, -i));
    stepsizes.push_back(5.0 * pow(10, -i));
    stepsizes.push_back(10.0 * pow(10, -i));
  }

  // begin fid sims
  for (auto &step : stepsizes) {

    dt = step;
    sim::dt_integration = dt;

    // Make FidFactory
    FidFactory ff;
    ff.SimulateFid(wf, tm);

    pt.Fill();

  } // grad

  pf.Write();
  pf.Close();

  return 0;
}
