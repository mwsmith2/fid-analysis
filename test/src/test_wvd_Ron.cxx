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
  int window = 200;
  double ti = 0.0;
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

  FidFactory ff;
  ff.SetLarmorFreq(ftruth + sim::mixdown_freq);
  // ff.IdealFid(wf_re, tm);

  for (int i = 0; i < fid_length; ++i) {
    wf_re.push_back(sin(45.101 * tm[i]*2*M_PI));
  }

  Fid myfid(wf_re, tm);
  
  std::vector<double> freq= dsp::fftfreq(fid_length, dt);
  //std::vector<double> wf_short = myfid.wf(); 

  // out.open("test_data/Ideal_fid_test.csv");
  // myfid.WriteMethodCsv(out);
  myfid.SavePlot("test_data/wvd_test_fid.pdf", "Test Waveform");
 
  // out.close();
  /*
  auto wf_im = dsp::hilbert(wf_re);
  arma::cx_vec wf(wf_re.size());

  for (int i = 0; i< wf_re.size(); ++i) {
    wf[i] = arma::cx_double(wf_re[i], wf_im[i]);
  }
  
  //wvd_cx is a Nx1 vector so in order to plot it I need to treat it like a matrix.
  std::vector<double> wvd (fid_length, 0.0);
  for (int i = 400; i< fid_length; i++) {
    arma::cx_vec wf_rc = dsp::acorrelation(wf, i, 200);
    arma::vec real_fft_acorr = arma::real(arma::fft(wf_rc));
    
    arma::uword max_bin;
    double max_val = real_fft_acorr.max(max_bin);

      // int max_bin = std::distance(real_fft_acorr.begin(),
      //                      std::max_element(real_fft_acorr.begin(), real_fft_acorr.end()));
    int i_range=0; int f_range = 0;
    int n = real_fft_acorr.size();
    for (int i =1 ; i<n+1; i++) {
      if (real_fft_acorr[max_bin+i]>.2*real_fft_acorr[max_bin]) {
        f_range= max_bin + i;
        n = f_range-i_range;
      }
      if ((max_bin-i>0)&(real_fft_acorr[max_bin-i]>.2*real_fft_acorr[max_bin])) {
        i_range = max_bin - i;
      }else if (i_range <0){
        i_range = 0;
        f_range = 2*max_bin;
      }
    }
    
    //create vector of normalized peak.
    double range = f_range - i_range;
    arma::vec peak = arma::zeros(range);
    double sum = 0;
    double dx = freq[1]-freq[0];
    for (int i=i_range; i<f_range+1; i++) {
      sum += .5*(real_fft_acorr[i+1]+real_fft_acorr[i])*dx;
    }
  
    for (int i = 0; i<peak.size(); i++) {
      peak[i]=real_fft_acorr[i_range+i]/sum;
    }
    
    //Now perform fit to Gaussian
    std::string gaussian("[2]*exp(-(x-[0])^2/(2*[1]^2))+[3]");
    TGraph gr3 = TGraph(range, &freq[i_range], &peak[0]);
    TF1 fit_func = TF1("fit_func", gaussian.c_str(), freq[i_range],freq[f_range-1]);
    
    fit_func.SetParameter(0, freq[max_bin]);
    fit_func.SetParameter(1, 2);
    fit_func.SetParameter(2, peak.max());
    fit_func.SetParameter(3, 0);
    
    gr3.Fit(&fit_func, "R");

    // std::cout<<"Center Frequency of FFT is: " << fit_func.GetParameter(0)/2<<endl;
    
    wvd.push_back(fit_func.GetParameter(0)/2);
    
  }
  */
  std::vector<double> wvd_f_t(fid_length, 0.0);
  
  wvd_f_t = dsp::wvd(wf_re , true, 200);
  
  std::cout<< "The centroid frequency is : "<< myfid.CalcCentroidFreq()<<std::endl;
  
  // auto real = arma::conv_to<std::vector<double>>::from(arma::real(res(arma::span(0,2500),2000)));
  // auto imag = arma::conv_to<std::vector<double>>::from(arma::imag(res(arma::span(0,2500), 2000)));
  // wvd_freq = dsp::wvd(wf_re, true, 200);

  //Plot wf(t)
  TMultiGraph mg;
  TCanvas c1;
  c1.Divide(2,1); 
  c1.cd(1);

  TGraph gr = TGraph(wf_re.size(), &tm[0], &wf_re[0]); gr.SetLineColor(kBlue);
  // TGraph gi = TGraph(wf.size()/2, &freq[0], &imag[0]); gi.SetLineColor(kRed);

  mg.Add(&gr, "cp"); gr.SetTitle("WVD frequencies");
  // mg.Add(&gi, "cp"); gi.SetTitle("Imaginary WVD component");

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
