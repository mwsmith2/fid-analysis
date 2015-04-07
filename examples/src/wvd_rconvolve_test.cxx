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
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  // declare variables
  int fid_length = 5000;
  double ti = -1.0;
  double dt = 0.001;
  double ftruth = 23.0;

  vector<double> wf_re;
  vector<double> tm;
  wf_re.reserve(fid_length);
  tm.reserve(fid_length);

  std::ofstream out;
  out.precision(10);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  for (int i = 0; i < fid_length; ++i) {
    wf_re.push_back(sin(40 * tm[i]));
  }

  auto wf_im = dsp::hilbert(wf_re);
  arma::cx_vec wf(wf_im.size());

  for (int i = 0; i < wf_re.size(); ++i) {
    wf[i] = arma::cx_double(wf_re[i], wf_im[i]);
  }
  

  auto wf_rc = dsp::rconvolve(wf, 4000);

  out.open("wvd_test_rc_real.txt");
  for (int i = 0; i < wf_rc.n_elem; ++i) {
    out << wf_rc[i].real() << ",";
  }
  out.close();

  out.open("wvd_test_rc_imag.txt");
  for (int i = 0; i < wf_rc.n_elem; ++i) {
    out << wf_rc[i].imag() << ",";
  }
  out.close();

  return 0;
}
