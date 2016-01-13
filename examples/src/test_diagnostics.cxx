/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This is a new test program for my Fid libraries 

\*===========================================================================*/


//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>

//--- other includes --------------------------------------------------------//
#include <armadillo>

//--- project includes ------------------------------------------------------//
#include "fid.h"

using std::vector;
using std::cout;
using std::endl;
using namespace fid;

int main(int argc, char** argv)
{
  // set precision
  cout.precision(10);
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  std::ofstream out;
  out.precision(10);

  // declare variables
  int fid_length = sim::num_samples;
  double ti = -1.0;
  double dt = 0.001;
  double ftruth = 23.0;

  vector<double> wf(0.0, fid_length);
  vector<double> tm(0.0, fid_length);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  FidFactory ff;
  ff.SetFreqLarmor(ftruth + sim::mixdown_freq);
  ff.IdealFid(wf, tm);

  Fid my_fid(wf, tm);

  my_fid.PrintDiagnosticData();
  my_fid.DumpDiagnosticData("test/");

  return 0;
}
