/*===========================================================================*\

Author: Ronaldo Ortez
Email:  supron00@gmail.com
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
#include "TGraph.h"
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
  int window = 60;
  double ti = -1.0;
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
  
  //  FidFactory ff;
  // ff.SetLarmorFreq(ftruth + sim::mixdown_freq);
  //ff.IdealFid(wf_re, tm);

  for (int i = 0; i < fid_length; i++) {
    if (tm[i]<0 ){ 
      wf_re.push_back(0.0);
    } else {
      wf_re.push_back(cos(45.101 * tm[i]*2*M_PI));
    }
  }

  Fid my_fid(wf_re, tm);

  // Extract Sinusoidal "FID" portion of the signal.
  unsigned int fid_begin = my_fid.i_wf();
  int N = wf_re.size();
  cout<< "fid begins at " <<fid_begin<<endl;

  // Produce FID centered around zero.
  std::vector<double> wf_cen(N, 0.0);
  cdouble r = cdouble(1.0, 0.0);
  std::vector<double> wf_pi = dsp::hilbert_rot(wf_re,r);
  std::transform( wf_re.begin(), wf_re.end(), wf_pi.begin(), wf_cen.begin(),
                  [](double re, double cen) {return re - std::abs(re-cen);});

  // Make armadillo analytic version of wf_re
  auto wf_im = dsp::hilbert(wf_cen);
  arma::cx_vec wf(wf_im.size());

  for (int i = 0; i < wf_cen.size(); i++) {
    wf[i] = arma::cx_double(wf_cen[i], wf_im[i]);
  }

  //Compute auto-correlation function.
  arma::cx_vec wf_rc = dsp::acorrelation(wf, fid_begin+10 ,window);

  //Print out text files of auto-correlation function
  out.open("test_data/acorr_test_real.txt");
  for (int i = 0; i < N; ++i) {
    out << wf_cen[i] << ",";
  }
  out.close(); 

  out.open("test_data/acorr_test_imag.txt");
  for (int i = 0; i < N; ++i) {
    out << wf_im[i] << ",";
  }
  out.close();
  cout<<"made text files of modified fid"<<endl;

  // Compute the power spectrum of the waveform.
  auto fft_wf = dsp::psd(wf_re);
  auto acorr = arma::conv_to<std::vector<cdouble>>::from(wf_rc);
  auto fft_acorr = dsp::fft(acorr);//returns vector of size N/2

  // Define vectors to draw on first part of canvas.
  auto real = arma::conv_to<std::vector<double>>::from(arma::real(wf));
  auto imag = arma::conv_to<std::vector<double>>::from(arma::imag(wf));

  TMultiGraph mg;
  TGraph gr = TGraph(100, &tm[fid_begin], &real[fid_begin]);gr.SetLineColor(kBlue);
  TGraph gi = TGraph(100, &tm[fid_begin], &imag[fid_begin]);gi.SetLineColor(kRed);
  TCanvas c1; 
  c1.Divide(2,2);
  c1.cd(1);  
  mg.Add(&gr, "cp");
  mg.Add(&gi, "cp"); 

  mg.Draw("ap");

  cout<< "made it past the first part of the canvas"<< endl;

  //Define vectors to draw on the second part of the canvas.
  c1.cd(2);
  auto fft_acf = arma::fft(wf_rc);
  std::vector<double> real_fft = arma::conv_to<std::vector<double>>::from(arma::real(fft_acf));
  std::vector<double> imag_fft = arma::conv_to<std::vector<double>>::from(arma::imag(fft_acf));
  std::vector<double> freqs = dsp::fftfreq(wf.size(), dt);

  TMultiGraph mg1;
  TGraph gr1 = TGraph(N/2, &freqs[0], &real_fft[0]);gr1.SetLineColor(kBlue);
  TGraph gi1 = TGraph(N/2, &freqs[0], &imag_fft[0]);gi1.SetLineColor(kRed);

  mg1.Add(&gr1, "cp");
  mg1.Add(&gi1, "cp");

  // TPad *p1 = (TPad *)(c1.cd(2));
  // p1->SetLogy();

  mg1.Draw("a");
  
  cout<<"made it past the second part of the cavas." << endl;

  //Define vectors to draw on the third part of the canvas.
  std::vector<double> fft_real;// (fid_length/2, 0.0);
  fft_real.reserve(fid_length/2);
  std::vector<double> fft_imag;//
  fft_imag.reserve(fid_length/2);
  for (int i =0;i < fid_length/2; i++) {
    fft_real.push_back(std::real(fft_wf[i]));
    fft_imag.push_back(std::imag(fft_wf[i]));
  }

  TMultiGraph mg2;
  TGraph gr2 = TGraph(N/50, &tm[0], &wf_cen[0]);gr2.SetLineColor(kBlue);
  TGraph gi2 = TGraph(N/50, &tm[0], &wf_im[0]);gi2.SetLineColor(kRed);

  mg2.Add(&gr2, "cp");
  mg2.Add(&gi2, "cp");

  // TPad *p2 = (TPad *)(c1.cd(3));
  // p2->SetLogy();

  c1.cd(3);
  mg2.Draw("ap");
 
  cout<< "made it past the third part of the canvas"<<endl;

  //Now find centroid by Gaussian fit
  int max_bin = std::distance(real_fft.begin(),
                              std::max_element(real_fft.begin(), real_fft.end()));
  int i_range=0; int f_range = 0;
  int n = real_fft.size();
  std::cout<<"peak bin is: " << max_bin<<std::endl;

  for (int i =1 ; i<n+1; i++) {
    if (real_fft[max_bin+i]>.1*real_fft[max_bin]) {
      f_range= max_bin + i+1;
      n = f_range-i_range;
    }
    if ((real_fft[max_bin-i]>.1*real_fft[max_bin])&(max_bin-i>0)) {
      i_range = max_bin - i;
      //   n=f_range-i_range;
    }else if (i_range <0){
      i_range = 0;
      f_range = 2*max_bin;
    }
    //  n=f_range-i_range;
  }

  //Remove offset for better fit.
  double range = f_range - i_range;
  arma::vec peak = arma::zeros(range);
  double sum = 0;
  double dx = freqs[1]-freqs[0];
  for (int i=i_range; i<f_range+1; i++) {
    sum += .5*(real_fft[i+1]+real_fft[i])*dx;
  }
  
  for (int i = 0; i<peak.size(); i++) {
    peak[i]=real_fft[i_range+i]/sum;
  }
 
  TMultiGraph mg3;
  TGraph gr3 = TGraph(range, &freqs[i_range], &peak[0]);

  mg3.Add(&gr3 , "cp");

  c1.cd(4);
  mg3.Draw("ap");
  
  std::cout<<" i_range is:  " << i_range<< "  f_range is:  "<<f_range<< std::endl;

  std::string gaussian("[2]*exp(-(x-[0])^2/(2*[1]^2))+[3]");
  TF1 fit_func = TF1("fit_func", gaussian.c_str(), freqs[i_range],freqs[f_range-1]);
  fit_func.SetParameter(0, freqs[max_bin]);
  fit_func.SetParameter(1, 2);
  fit_func.SetParameter(2, peak.max());
  fit_func.SetParameter(3, 0);

  gr3.Fit(&fit_func, "R");

  std::cout<<"Center Frequency of FFT is: " << fit_func.GetParameter(0)/2<<endl;
  std::cout<<"Center Frequency from Centroid is: "<<my_fid.CalcCentroidFreq();

  //Save canvas.
  c1.Print("test_data/acorr_test.pdf");

  return 0;
}
