/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   15/04/14

Detail: The program is meant to test the effects of linear field 
        gradients on the FID frequency extraction.  The sweep parameters
        are set in a separate config file here, but the user need not 
        rely on the config parameters.  All that needs to be done is 
        the defining of a gradient vector.

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
  // Load the configurable parameters
  if (argc > 1) {
    load_params(std::string(argv[1]));
  }

  // allocate some necessary parameters
  std::vector<double> wf;

  double final_time = sim::start_time + sim::num_samples*sim::delta_time;
  std::vector<double> tm = construct_range(sim::start_time, final_time, sim::delta_time);

  // csv output
  std::cout.precision(10);
  std::cout.setf(std::ios::fixed, std::ios::floatfield);

  // Make FidFactory
  FidFactory ff;
  ff.SimulateFid(wf, tm);

  FID my_fid(wf, tm);

  my_fid.SavePlot("test.png", "FID Test");

  fid::write_fid_file("test.fid", my_fid.tm(), my_fid.wf());

  return 0;
}
