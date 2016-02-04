/*===========================================================================*\

Author: Ronaldo Ortez
Email:  mwsmith2@uw.edu
Date:   1/23/16

Detail: This is a test program which test the auto-correlation function I wrote  

\*===========================================================================*/


//--- std includes ----------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
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
  cout.setf(std::ios::fixed, std:: ios::floatfield);

  // declare variables
  int fid_length = 5000;
  double ti = 0;
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
  
  //FidFactory::IdealFid(wf_re, tm);
  for (int i = 0; i < fid_length; ++i) {
    wf_re.push_back(sin(50 * tm[i]*6.28));
  }

  auto wf_im = dsp::hilbert(wf_re);
  arma::cx_vec wf(wf_im.size());

  for (int i = 0; i < wf_re.size(); ++i) {
    wf[i] = arma::cx_double(wf_re[i], wf_im[i]);
  }
  
  //Compute auto-correlation function.
  arma::cx_vec wf_rc = dsp::acorrelation(wf, 2000, 200);

  //Print out text files of auto-correlation function
  out.open("test_data/acorr_test_real.txt");
  for (int i = 0; i < wf_rc.n_elem; ++i) {
    out << wf_rc[i].real() << ",";
  }
  out.close(); 

  out.open("test_data/acorr_test_imag.txt");
  for (int i = 0; i < wf_rc.n_elem; ++i) {
    out << wf_rc[i].imag() << ",";
  }
  out.close();

  //Compute the power spectrum of the waveform.
  auto fft_wf = dsp::psd(wf_re);
  auto acorr = arma::conv_to<std::vector<cdouble>>::from(wf_rc);
  auto fft_acorr = dsp::fft(acorr);//returns vector of size N/2

  //Define vectors to draw on first part of canvas.
  auto real = arma::conv_to<std::vector<double>>::from(arma::real(wf_rc));
  auto imag = arma::conv_to<std::vector<double>>::from(arma::imag(wf_rc));

  TMultiGraph mg;
  TGraph gr = TGraph(wf.size(), &tm[0], &real[0]);gr.SetLineColor(kBlue);
  TGraph gi = TGraph(wf.size(), &tm[0], &imag[0]);gi.SetLineColor(kRed);
  TCanvas c1;
  c1.Divide(2,2);
  c1.cd(1);
  mg.Add(&gr, "cp");
  mg.Add(&gi, "cp");

  mg.Draw("ap");

  
  //Define vectors to draw on the second part of the canvas.
  c1.cd(2);
  auto fft_acf = arma::fft(wf_rc);
  std::vector<double> real_fft = arma::conv_to<std::vector<double>>::from(arma::real(fft_acf));
  std::vector<double> imag_fft = arma::conv_to<std::vector<double>>::from(arma::imag(fft_acf));
  std::vector<double> freqs = dsp::fftfreq(wf.size(), dt);

  TMultiGraph mg1;
  TGraph gr1 = TGraph(wf.size()/2, &freqs[0], &real_fft[0]);gr1.SetLineColor(kBlue);
  TGraph gi1 = TGraph(wf.size()/2, &freqs[0], &imag_fft[0]);gi1.SetLineColor(kRed);

  mg1.Add(&gr1, "cp");
  mg1.Add(&gi1, "cp");

  TPad *p1 = (TPad *)(c1.cd(2));
  p1->SetLogy();

  mg1.Draw("a");

  //Define vectors to draw on the third part of the canvas.
  std::vector<double> fft_real;// (fid_length/2, 0.0);
  fft_real.reserve(fid_length/2);
  std::vector<double> fft_imag;//
  fft_imag.reserve(fid_length/2);
  for (int i =0;i < fid_length/2; i++) {
    fft_real.push_back(std::real(fft_acorr[i]));
    fft_imag.push_back(std::imag(fft_acorr[i]));
  }

  TMultiGraph mg2;
  TGraph gr2 = TGraph(wf.size()/2, &freqs[0], &fft_real[0]);gr2.SetLineColor(kBlue);
  TGraph gi2 = TGraph(wf.size()/2, &freqs[0], &fft_imag[0]);gi2.SetLineColor(kRed);

  mg2.Add(&gr2, "cp");
  mg2.Add(&gi2, "cp");

  TPad *p2 = (TPad *)(c1.cd(3));
  p2->SetLogy();

  c1.cd(3);
  mg2.Draw("ap");

  //Save canvas.
  c1.Print("test_data/acorr_test.pdf");


  return 0;
}
