/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   2015/01/08

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
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  double grad_min = 0;
  double grad_max = 100;
  double dgrad = 1.0;
  char str[60];

  // some necessary parameters
  vec wf;
  vec tm;
  vec grads;
  vec grad_0;
  vec gradient;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  tm = construct_range(sim::start_time, final_time, sim::delta_time);

  grads = construct_range(grad_min, grad_max, dgrad);

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

