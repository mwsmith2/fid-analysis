/*===========================================================================*\

Author: Ronaldo Ortez
Email:  supron00@gmail.com
Date:   2/22/15

Detail: This is a new test program for the WVD technique. 

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
#include "TMultiGraph.h"
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
  int window = 80;
  double ti = 0.0;
  double dt = 0.001;
  double ftruth = 23.0;

  vector<double> wf_re;
  vector<double> tm;
  wf_re.reserve(fid_length);
  tm.reserve(fid_length);

  std::ofstream out;
  out.precision(10);

  // Create Sample waveform.
  for (int i = 0; i < fid_length; i++){
    tm.push_back(i * dt + ti);
  }

  FidFactory ff;
  ff.SetLarmorFreq(ftruth + sim::mixdown_freq);
  ff.IdealFid(wf_re, tm);

  /* for (int i = 0; i < fid_length; ++i) {
    wf_re.push_back(sin(45.101 * tm[i]*2*M_PI));
    }*/

  Fid myfid(wf_re, tm);

  // Extract FID portion of the signal.
  int fid_begin = myfid.i_wf();
  wf_re.erase(wf_re.begin(), wf_re.begin() + fid_begin);
  int N = wf_re.size();
  
  // Produce FID centered around zero.
  std::vector<double> wf_cen(N, 0.0);
  cdouble r = cdouble(1.0, 0.0);
  std::vector<double> wf_pi = dsp::hilbert_rot(wf_re,r);
  std::transform( wf_re.begin(), wf_re.end(), wf_pi.begin(), wf_cen.begin(),
                  [](double re, double cen) {return re - std::abs(re-cen);});

  // Get frequency dimensionality and save test waveform
  std::vector<double> freq= dsp::fftfreq(fid_length, dt);
  myfid.SavePlot("test_data/wvd_test_fid.pdf", "Test Waveform");
 

  // Extract frequency with time using wvd technique
  std::vector<double> wvd_f_t(fid_length, 0.0);
  wvd_f_t = dsp::WvdFreqExt(wf_cen, true, window);
  
  std::cout<< "The centroid frequency is : "<< myfid.CalcCentroidFreq()<<std::endl;
  

  //Plot wf(t)
  TMultiGraph mg;
  TCanvas c1;
  c1.Divide(2,1); 
  c1.cd(1);

  TGraph gr = TGraph(wf_re.size()/50, &tm[0], &wf_cen[0]); gr.SetLineColor(kBlue);
  TGraph gi = TGraph(wf_re.size()/50, &tm[0], &wf_re[0]); gi.SetLineColor(kRed);

  mg.Add(&gr, "cp"); gr.SetTitle("WVD frequencies");
  mg.Add(&gi, "cp"); gi.SetTitle("Imaginary WVD component");

  mg.Draw("al");


  //Plot WVD(f,t)
  c1.cd(2);
  TMultiGraph mg1;
  TGraph gr1 = TGraph(wvd_f_t.size()-2*window, &tm[window], &wvd_f_t[window]);
  mg1.Add(&gr1 , "cp");

  mg1.Draw("al");
 
 
  c1.Print("test_data/wvd_test.pdf");
  
  return 0;
}
