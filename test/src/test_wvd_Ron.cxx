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
  double ti = 0.0;
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

  FidFactory ff;
  ff.SetLarmorFreq(ftruth + sim::mixdown_freq);
  ff.IdealFid(wf, tm);

  Fid myfid(wf, tm);
  
  std::vector<double> freq= dsp::fftfreq(fid_length, dt);

  out.open("test_data/Ideal_fid_test.csv");
  myfid.WriteMethodCsv(out);
  myfid.SavePlot("test_data/wvd_test_fid.pdf", "Test Waveform");
  // for (int i = 0; i < wf.size(); ++i) {
  // wf[i] = sin(40 * tm[i]);
  // }

  /* out.open("wvd_test_wf.txt");
  for (auto it = wf.begin(); it != wf.end(); ++it) {
    out << *it << ", ";
  }
  out.close();

  auto wf_im = dsp::hilbert(wf);

  out.open("wvd_test_wf_im.txt");
  for (auto it = wf_im.begin(); it != wf_im.end(); ++it) {
    out << *it << ", ";
    }*/
  out.close();
  
  //wvd_cx is a NxN matrix so in order to plot it I need to treat it like a matrix.
  auto res = dsp::wvd_cx(wf, true, 200);
  // res.quiet_save("wvd_test.txt", arma::csv_ascii);
  
  std::cout<< "The centroid frequency is : "<< myfid.CalcCentroidFreq()<<std::endl;
  
  auto real = arma::conv_to<std::vector<double>>::from(arma::real(res(arma::span(0,2500),2000)));
  auto imag = arma::conv_to<std::vector<double>>::from(arma::imag(res(arma::span(0,2500), 2000)));

  // Fid my_wvd_real(subReal,freq);
  // Fid my_wvd_imag(imag,freq);
  
  TMultiGraph mg;
  TCanvas c1;

  TGraph gr = TGraph(wf.size()/2, &freq[0], &real[0]);
  gr.SetLineColor(kBlue);
  TGraph gi = TGraph(wf.size()/2, &freq[0], &imag[0]); gi.SetLineColor(kRed);

  // gr = TGraph(wf.size()/2, &freq[0], &subReal[0]);
  // gi = TGraph(wf.size()/2, &freq[0], &imag[0]);

  mg.Add(&gr, "cp"); gr.SetTitle("Real WVD component");
  mg.Add(&gi, "cp"); gi.SetTitle("Imaginary WVD component");

  mg.Draw("a");
  // gr.GetXaxis();
  //gr.SetLimits(0.0,2.5);

  c1.Print("test_data/wvd_test.pdf");
  
  /* gi.Draw();
  c1.Update();
  c1.Print("test_data/wvd_imag_test.pdf");*/
    // my_wvd_real.SavePlot("test_data/wvd_real_test.pdf" , "Real WVD component");
    // my_wvd_imag.SavePlot("test_data/wvd_imag_test.pdf", "Image WVD compoent");

  return 0;
}
