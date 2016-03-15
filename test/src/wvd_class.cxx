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
  double ti = -0.5;
  double dt = 0.001;
  double ftruth = 100.03;

  bool upsample = false;
  int window = 125;
  int N= fid_length;

  vector<double> wf;
  vector<double> tm;
  wf.reserve(fid_length);
  tm.reserve(fid_length);

  std::ofstream out;
  out.precision(10);

  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }
  wf.resize(N, 0.0);
  //  FidFactory ff;
  //ff.SetLarmorFreq(ftruth + sim::mixdown_freq);
  //ff.IdealFid(wf, tm);

  for (int i = 0; i < wf.size(); ++i) {
    if (tm[i]>0) {
      wf[i] = sin((40.1 * tm[i]+ .1*.2*sin(tm[i]/.2))*2*M_PI);// sinusoid with sinusoidal
                                                              // frequency variation
    }
  }
  
  WvdFid mywvd(wf, tm, true, window); 

  //Plot results
  vector<double> wvd_f = mywvd.wvd_f();
  int M = wvd_f.size();

  uint fid_start = mywvd.i_wf();
  
  TMultiGraph mg; 
  TGraph gr = TGraph(M-2*fid_start, &tm[fid_start+10], &wvd_f[fid_start+10]);gr.SetLineColor(kBlue);
  TCanvas c1;
  mg.Add(&gr, "cp");
  mg.SetTitle("WVD Freq Extraction; time[ms]; freq[kHz]");
  mg.Draw("ap");
  
  c1.Print("test_data/wvd_class_freq_ext.pdf");
  // mywvd.SaveWvd("test_data/wvd_class_freq_ext.pdf", "WVD Freq Extraction");

  return 0;
}
