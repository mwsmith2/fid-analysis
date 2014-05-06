/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   15/04/14

Detail: The program is meant to generate a set of example gradient FIDs in
        plaintext format.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>
#include <iostream>        

//--- other includes --------------------------------------------------------//
#include "TFile.h"
#include "TTree.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_class.h"
#include "fid_sim.h"
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
  vec grads;
  vec grad_0;
  vec gradient;

  char str[60];

  construct_time_vector(len_fids, i_time, d_time, tm);
  grads = construct_sweep_range(grad_range);

  // Make FidFactory
  GradientFidFactory gff;
  construct_quadratic_gradient(20, grad_0);

  // csv output
  std::ofstream out;
  out.precision(10);
  out.setf(std::ios::fixed, std::ios::floatfield);

  // begin sweeps
  for (auto g: grads){

    sprintf(str, "data/test_fids/quad_%03dppb.fid", (int)g);
    out.open(str);

    gradient.resize(0);

    for (auto val : grad_0){
      gradient.push_back(val * g);
    }

    gff.ConstructFid(gradient, wf);

    for (int i = 0; i < tm.size(); ++i){

      out << tm[i] << " " << wf[i] << endl;

    }

    out.close();

  } // grad

  return 0;
}

