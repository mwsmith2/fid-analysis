/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This is a new test program for my FID libraries 

\*===========================================================================*/


//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;


//--- project includes ------------------------------------------------------//
#include "fid.h"
#include "fid_class.h"


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

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  fid::getIdealFID(wf, tm, ftruth);

  fid::FID my_fid(wf, tm);

  cout << my_fid.CalcZeroCountFreq() << endl;
  cout << my_fid.CalcCentroidFreq() << endl;
  cout << my_fid.CalcAnalyticalFreq() << endl;
  cout << my_fid.CalcLorentzianFreq() << endl;
  cout << my_fid.CalcSoftLorentzianFreq() << endl;
  cout << my_fid.CalcExponentialFreq() << endl;
  cout << my_fid.CalcPhaseFreq() << endl;
  cout << my_fid.CalcSinusoidFreq() << endl;

  return 0;
}