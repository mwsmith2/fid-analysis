/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   15/04/14

Detail: The program is meant to test the effects of linear field 
        gradients on the Fid frequency extraction.  The sweep parameters
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
  std::vector<double> tm = time_vector(); // default params in ::fid

  // csv output
  std::cout.precision(10);
  std::cout.setf(std::ios::fixed, std::ios::floatfield);

  // Make FidFactory
  FidFactory ff;
  ff.SimulateFid(wf, tm);

  Fid my_fid(wf, tm);

  my_fid.SavePlot("test.png", "Fid Test");
  my_fid.SaveData("test.fid");

  return 0;
}
