/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   2015/01/08

Detail: The program is meant to generate a set of example gradient FIDs in
        plaintext format.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <iostream>        
#include <fstream>
#include <vector>
#include <string>
using std::cout;
using std::endl;        

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
  std::vector<double> wf;
  std::vector<double> tm;
  std::vector<double> grads;
  std::vector<double> grad_0;
  std::vector<double> gradient;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  tm = construct_range(sim::start_time, final_time, sim::delta_time);

  grads = construct_range(grad_min, grad_max, dgrad);

  // Make FidFactory
  FidFactory ff;
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

    ff.GradientFid(gradient, wf);

    for (int i = 0; i < tm.size(); ++i){

      out << tm[i] << " " << wf[i] << endl;

    }

    out.close();

  } // grad

  return 0;
}

