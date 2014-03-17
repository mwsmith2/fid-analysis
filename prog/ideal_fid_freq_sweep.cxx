/*===========================================================================*\

Author: Matthias W. Smith
Email:  mwsmith2@uw.edu
Date:   11/02/14

Detail: This program

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
  out.open("ideal_freq_sweep_data.csv");

  // set precision
  cout.precision(10);
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  // declare variables
  int fid_length = 10000;
  int nfids = 1000;
  double ti = -1.0;
  double dt = 0.01;

  double fi = 23.0;
  double ff = 23.2;
  double df = 0.001;

  vector<double> wf;
  vector<double> tm;
  wf.reserve(fid_length);
  tm.reserve(fid_length);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  for (double f = fi; f <= ff; f += df){

    cout << "Running for frequency " << f << ".\n";

    for (int i = 0; i < nfids; i++){

      fid::getIdealFID(wf, tm, f);

      fid::FID my_fid(wf, tm);

      out << f << ", ";
      out << my_fid.CalcZeroCountFreq() << ", ";
      out << 0.0 << ", ";
      out << my_fid.CalcCentroidFreq() << ", ";
      out << 0.0 << ", ";
      out << my_fid.CalcAnalyticalFreq() << ", ";
      out << my_fid.chi2() << ", ";
      out << my_fid.CalcLorentzianFreq() << ", ";
      out << my_fid.chi2() << ", ";
      out << my_fid.CalcSoftLorentzianFreq() << ", ";
      out << my_fid.chi2() << ", ";
      out << my_fid.CalcExponentialFreq() << ", ";
      out << my_fid.chi2() << ", ";
      out << my_fid.CalcPhaseFreq() << ", ";
      out << my_fid.chi2() << ", ";
      out << my_fid.CalcSinusoidFreq() << ", ";
      out << my_fid.chi2() << endl;

    }
  }

  out.close();
  return 0;
}