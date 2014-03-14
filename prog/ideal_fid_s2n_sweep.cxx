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
  std::ofstream out;
  out.precision(10);
  out.setf(std::ios::fixed, std:: ios::floatfield);
  out.open("ideal_s2n_sweep_data.csv");

  // set precision
  cout.precision(10);
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  // declare variables
  int fid_length = 10000;
  int nfids = 100;
  double ti = -1.0;
  double dt = 0.001;

  double fi = 23.0;
  double ff = 23.2;
  double df = 0.01;

  double phi = 0.0;
  double s2n_i = 50.0;
  double s2n_f = 200.0;
  double d_s2n = 10.0;

  vector<double> wf;
  vector<double> tm;
  wf.reserve(fid_length);
  tm.reserve(fid_length);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  for (double f = fi; f <= ff; f += df){

    cout << "Running for frequency " << f << ".\n";

    for (double s2n = s2n_i; s2n <= s2n_f; s2n += d_s2n){

      for (int i = 0; i < nfids; i++){

          fid::getIdealFID(wf, tm, f, phi, s2n);

          fid::FID my_fid(wf, tm);

          out << f << ", ";
          out << s2n << ", ";
          out << my_fid.CalcZeroCountFreq() << ", ";
          out << my_fid.CalcCentroidFreq() << ", ";
          out << my_fid.CalcAnalyticalFreq() << ", ";
          out << my_fid.CalcLorentzianFreq() << ", ";
          out << my_fid.CalcSoftLorentzianFreq() << ", ";
          out << my_fid.CalcExponentialFreq() << ", ";
          out << my_fid.CalcPhaseFreq() << ", ";
          out << my_fid.CalcSinusoidFreq() << endl;
      } // nfids
    } // phi
  } // tm

  out.close();
  return 0;
}