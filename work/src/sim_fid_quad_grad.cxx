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

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_class.h"
#include "fid_utilities.h"
#include "fid_sim.h"

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

  fid::ConstructTimeVector(len_fids, i_time, d_time, tm);
  grads = ConstructSweepRange(grad_range);

  // Make FidFactory
  GradientFidFactory gff;
  ConstructQuadraticGradient(20, grad_0);

  // csv output
  std::ofstream out;
  out.precision(10);
  out.setf(std::ios::fixed, std::ios::floatfield);

  // begin sweeps
  for (auto g: grads){

    cout << "gradient strength " << g << endl;

    for (int i = 0; i < num_fids; ++i){

      gradient.resize(0);

      for (auto val : grad_0){
        gradient.push_back(val * g);
      }

      gff.ConstructFid(gradient, wf);
      FID my_fid(wf, tm);

      CalcFreqSaveCsv(my_fid, out);

      if (i == 0){
        static char str[60];
        sprintf(str, "data/fig/fid_quad_grad_%03dppb.pdf", (int)g);
        DrawFID(my_fid, str, string("Test FID")); 
      }
    }
  } // grad

  return 0;
}
