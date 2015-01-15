/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   15/04/14

Detail: The program is meant to test the effects of field gradients
        on the FID frequency extraction.  The sweep parameters are
        set in a separate config file here, but the user need not rely
        on the config parameters.  All that needs to be done is the
        defining of a gradient vector.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <fstream>

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
  vec wf;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  vec tm = construct_range(sim::start_time, final_time, sim::delta_time);

  vec grads = construct_range(grad_min, grad_max, dgrad);
  vec grad_0;
  construct_quadratic_gradient(20, grad_0);

  // Make FidFactory
  FidFactory ff;

  // csv output
  std::ofstream out;
  out.precision(10);
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.open("data/sim_fid_quad_grad_data.csv");

  // begin sweeps
  for (auto g: grads){

    cout << "gradient strength " << g << endl;

    for (int i = 0; i < num_fids; ++i){

      vec gradient;

      for (auto val : grad_0){
        gradient.push_back(val * g);
      }

      ff.GradientFid(gradient, wf);
      FID my_fid(wf, tm);

      calc_freq_write_csv(my_fid, out);

      if (i == 0){
        static char str[60];
        sprintf(str, "data/fig/fid_quad_grad_%03dppb.pdf", (int)g);
        draw_fid(my_fid, str, string("Test FID")); 
      }
    }
  } // grad

  return 0;
}

