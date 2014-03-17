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
  double ftruth = 23.456789;
  int nfids = 10000;

  vector<double> wf;
  vector<double> tm;
  wf.reserve(fid_length);
  tm.reserve(fid_length);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  std::ofstream out;
  out.open("ideal_fid_test_data.csv");
  out.precision(10);
  out.setf(std::ios::fixed, std:: ios::floatfield);

  for (int i = 0; i < nfids; i++){

    fid::getIdealFID(wf, tm, ftruth);
    fid::FID my_fid(wf, tm);

    out << ftruth << ", ";
    out << my_fid.CalcZeroCountFreq() << ", ";
    out << my_fid.CalcCentroidFreq() << ", ";
    out << my_fid.CalcAnalyticalFreq() << ", ";
    out << my_fid.CalcLorentzianFreq() << ", ";
    out << my_fid.CalcSoftLorentzianFreq() << ", ";
    out << my_fid.CalcExponentialFreq() << ", ";
    out << my_fid.CalcPhaseFreq() << ", ";
    out << my_fid.CalcSinusoidFreq() << endl;
  }
  
  out.close();
  return 0;
}