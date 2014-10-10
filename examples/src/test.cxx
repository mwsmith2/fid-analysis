/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This is a new test program for my FID libraries 

\*===========================================================================*/


//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;

//--- other includes --------------------------------------------------------//
#include <armadillo>

//--- project includes ------------------------------------------------------//
#include "fid.h"
using namespace fid;


int main(int argc, char** argv)
{
  // set precision
  cout.precision(10);
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  // declare variables
  int fid_length = 10000;
  double ti = -1.0;
  double dt = 0.001;
  double ftruth = 23.0;

  vector<double> wf;
  vector<double> tm;
  wf.reserve(fid_length);
  tm.reserve(fid_length);

  std::ofstream out;
  out.precision(10);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  ideal_fid(wf, tm, ftruth);

  FID my_fid(wf, tm);

  auto res = dsp::wvd(wf);
  res.quiet_save("wvd_test.dat", arma::csv_ascii);

  return 0;
}
