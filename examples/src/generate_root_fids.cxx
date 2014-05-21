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
#include "TString.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_class.h"
#include "fid_utilities.h"
#include "fid_sim.h"

using namespace fid;
using namespace fid::sim;

int main(int argc, char **argv)
{
  // initialize the configurable parameters
  load_params(argc, argv);
  len_fids = num_points / reduction; // from fid::sim

  // some necessary parameters
  vec wf;
  vec tm;
  wf.reserve(len_fids);
  tm.reserve(len_fids);
  double delta_b;
  double max_grad = 500; // in ppm
  double freq_0 = sim::freq_larmor;

  construct_time_vector(len_fids, i_time, d_time, tm);

  // Make FidFactory
  FidFactory ff;

  // Set up the ROOT tree to hold the results
  TFile pf("test.root", "recreate");
  TTree pt("t", "FID Tree");

  pt.Branch("db", &delta_b, "delta_b/D");
  pt.Branch("fid", &wf[0], TString::Format("fid[%d]/D", len_fids));
  pt.Branch("tm", &tm[0], TString::Format("time[%d]/D", len_fids));

  // begin grad sims
  for (int i = -1 * num_fids / 2; i < num_fids / 2 + 1; ++i){

    delta_b = max_grad * 2 * i / num_fids;

    sim::freq_larmor = 1.0e-6 * delta_b + freq_0;

    ff.SimulateFid(wf, tm);

    pt.Fill();

  } // grad

  pf.Write();
  pf.Close();

  return 0;
}
