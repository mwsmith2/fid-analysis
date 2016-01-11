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
#include <cmath>
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
  cout.setf(std::ios::fixed, std::ios::floatfield);

  // declare variables
  int fid_length = 5000;
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

  sim::freq_larmor = ftruth + sim::freq_ref;

  FidFactory ff;
  ff.IdealFid(wf, tm);
  
  for (int i = 0; i < wf.size(); ++i) {
    wf[i] = sin(40 * tm[i]);
  }

  out.open("wvd_test_wf.txt");
  for (auto it = wf.begin(); it != wf.end(); ++it) {
    out << *it << ", ";
  }
  out.close();

  auto wf_im = dsp::hilbert(wf);

  out.open("wvd_test_wf_im.txt");
  for (auto it = wf_im.begin(); it != wf_im.end(); ++it) {
    out << *it << ", ";
  }
  out.close();

  auto res = dsp::wvd(wf);
  res.quiet_save("wvd_test.txt", arma::csv_ascii);

  return 0;
}
