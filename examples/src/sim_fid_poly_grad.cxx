/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   2015/01/08

Detail: The program is meant to test the effects of field gradients
        on the FID frequency extraction.  The sweep parameters are
        set in a separate config file here, but the user need not rely
        on the config parameters.  All that needs to be done is the
        defining of a gradient vector.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
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
  char str[60];
  int npoints = 21;
  int num_fids = 100;

  // some necessary parameters
  std::vector<double> wf;
  std::vector<double> tm;
  std::vector<double> grads;
  std::vector<double> grad_0;
  std::vector<double> gradient;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  tm = construct_range(sim::start_time, final_time, sim::delta_time);

  // Make FidFactory
  FidFactory ff;

  // construct a normalized, centered polynomial gradient
  sprintf(str, "pol%d", grad::poln_order);
  TF1 f1("f1", str);

  for (auto &val : grad::poln_coefs){
    static int count = 0;
    f1.SetParameter(count++, val);
  }

  // Now set the values in the gradient
  double dx = 10.0 / (npoints - 1);
  for (double x = -5.0; x <= 5.0; x += 1.0){
    grad_0.push_back(f1.Eval(x));
  }

  // now center and normalize vector
  double avg = 
    std::accumulate(grad_0.begin(), grad_0.end(), 0.0) / grad_0.size();
  for (int i = 0; i < grad_0.size(); i++){
    grad_0[i] -= avg;
  }

  // normalize by largest value
  double max = *std::max_element(grad_0.begin(), grad_0.end());
  for (int i = 0; i < grad_0.size(); i++){
    grad_0[i] /= max;
  }

  // csv output
  std::ofstream out;
  out.precision(10);
  out.setf(std::ios::fixed, std::ios::floatfield);
  sprintf(str, "data/sim_fid_pol%d_grad_data.csv", grad::poln_order);
  out.open(str);

  // begin sweeps
  for (auto g: grads){

    cout << "gradient strength " << g << endl;

    for (int i = 0; i < num_fids; ++i){

      gradient.resize(0);

      for (auto val : grad_0){
        gradient.push_back(val * g);
      }

      ff.GradientFid(gradient, wf);
      FID my_fid(wf, tm);

      calc_freq_write_csv(my_fid, out);

      if (i == 0){
        sprintf(str, "data/fig/fid_pol%d_grad_%03dppb.pdf", grad::poln_order, (int)g);
        draw_fid(my_fid, str, std::string("Test FID")); 
      }
    }
  } // grad

  return 0;
}

