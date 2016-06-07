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
  int fid_length = 10000;
  double ti = -0.5;
  double dt = 0.001;
  double ftruth = 101.030;

  bool upsample = true;
  int window = 500;
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
  // FidFactory ff;
  // ff.SetLarmorFreq(ftruth + sim::mixdown_freq);
  // ff.IdealFid(wf, tm);

  for (int i = 0; i < wf.size(); ++i) {
    if (tm[i]>0) {
      wf[i] = sin(43.101 * tm[i]*2*M_PI+ .1*sin(tm[i]*2*M_PI));// sinusoid with sinusoidal
      // frequency variation
    }
  }
  
  int fid_stop = 1500;
  WvdFid mywvd(wf, tm); 
  double wvd_f0 = mywvd.WvdZeroFreq(fid_stop);

  //Plot results
  vector<double> wvd_f = mywvd.wvd_f();
  int M = wvd_f.size();
  vector<cdouble> analytic_fid = mywvd.AnalyticWf();
  cout << "The quoted size of analytic_wf_ is: " << analytic_fid.size()<<endl;
  vector<double> re_fid(analytic_fid.size(), 0.0);
  vector<double> im_fid(analytic_fid.size(), 0.0);
  for (int i = 0; i< analytic_fid.size();i++) {
    re_fid[i] = real(analytic_fid[i]);
    im_fid[i] = imag(analytic_fid[i]);
  }

  uint fid_start = mywvd.i_wf();
  
  cout<< " The size of the output frequency extraction vector is: "<<M
      << " and beginning of the fid is: "<<fid_start<<endl;

  // Plot things
  TMultiGraph mg; 
  TGraph gr = TGraph(M-400, &tm[fid_start+200], &wvd_f[200]);gr.SetLineColor(kBlue);
  TCanvas c1;
  c1.Divide(2,2);
  c1.cd(1);
  mg.Add(&gr, "P");
  mg.SetTitle("WVD Freq Extraction; time[ms]; freq[kHz]");
  mg.Draw("ap");

  //Second pad of Canvas
  TMultiGraph mg1;
  TGraph wf_re = TGraph(re_fid.size()/50, &tm[fid_start], &re_fid[0]);wf_re.SetLineColor(kBlue);
  TGraph wf_im = TGraph(im_fid.size()/50, &tm[fid_start], &im_fid[0]);wf_im.SetLineColor(kRed);
  
  c1.cd(2);
  mg1.Add(&wf_re, "cp");
  mg1.Add(&wf_im, "cp");
  mg1.Draw("ap");

  // Third pad of Canvas
  auto extrema = std::minmax_element(wvd_f.begin()+100, wvd_f.begin()+fid_stop);
  TH1D h1 = TH1D("h1", "Extracted Frequency Distribution", 
                 1000,*extrema.first, *extrema.second);
  for (int i = 100; i<fid_stop; i++) {
    h1.Fill(wvd_f[i]);
  }
  
  c1.cd(3);
  h1.Draw();
  
  cout<<" the WVD zero time frequency is: "<< wvd_f0<<endl;

  c1.Print("test_data/wvd_class_freq_ext.pdf");
  // mywvd.SaveWvd("test_data/wvd_class_freq_ext.pdf", "WVD Freq Extraction");

  return 0;
}
