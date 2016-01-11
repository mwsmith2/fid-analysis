/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   15/04/14

Detail: The program is meant to test the effects of linear field gradients on the Fid frequency extraction.  The sweep parameters are set in a separate config file here, but the user need not rely on the config parameters.  All that needs to be done is the defining of a gradient vector.

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
  // initialize the configurable parameters
  int num_fids = 100;

  double grad_min = 0;
  double grad_max = 100;
  double dgrad = 1.0;

  // allocate some necessary parameters
  std::vector<double> wf;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  std::vector<double> tm = construct_range(sim::start_time, final_time, sim::delta_time);

  std::vector<double> grads = construct_range(grad_min, grad_max, dgrad);
  std::vector<double> grad_0;
  grad_0 = normalized_gradient(20, 1);

  // Make FidFactory
  FidFactory ff;

  // csv output
  std::ofstream out;
  out.precision(10);
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.open("data/sim_fid_lin_grad_data.csv");

  // begin sweeps
  for (auto g: grads){

    cout << "gradient strength " << g << endl;

    for (int i = 0; i < num_fids; ++i){

      std::vector<double> gradient;

      for (auto val : grad_0){
        gradient.push_back(val * g);
      }

      ff.GradientFid(gradient, wf);
      Fid my_fid(wf, tm);

      my_fid.WriteFreqCsv(out);

      if (i == 0){
        static char str[60];
        sprintf(str, "data/fig/fid_lin_grad_%03dppb.pdf", (int)g);
        my_fid.SavePlot(str); 
      }
    }
  } // grad

  return 0;
}
