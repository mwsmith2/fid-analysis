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
  int fid_length = 10000;
  int window =50;
  double ti = -1.0;
  double dt = 0.001;
  double ftruth = 43.0; 

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
  // ff.IdealFid(wf_re, tm);

  /*  for (int i = 0; i < fid_length; i++) {
    if (tm[i]<0 ){ 
      wf_re.push_back(0.0);
    } else {
      wf_re.push_back(cos(43.101 * tm[i]*2*M_PI));
    }
  }

  Fid my_fid(wf_re, tm);*/

  Fid my_fid("../fids/fid_ch00_00001.fid");
  wf_re = my_fid.wf();
  tm = my_fid.tm();
            
  // Extract Sinusoidal "FID" portion of the signal.
  unsigned int fid_begin = my_fid.i_wf();
  // fid_begin = 1000;
  wf_re.erase(wf_re.begin(), wf_re.begin() + fid_begin);
  int N = wf_re.size();
  cout<< "fid begins at " <<fid_begin<<endl;

  vector<double> wf_fil = dsp::hilbert_rot(wf_re,(0.0,0.0),10);
  vector<double> wf_fil_rot = dsp::hilbert_rot(wf_fil, cdouble(0.0, 1.0));

  // Produce FID centered around zero.
  std::vector<double> wf_cen(N, 0.0);
  /* cdouble r = cdouble(1.0, 0.0);
  std::vector<double> wf_pi = dsp::hilbert_rot(wf_re,r);
  std::transform( wf_re.begin(), wf_re.end(), wf_pi.begin(), wf_cen.begin(),
                  [](double re, double cen) {return re - std::abs(re-cen);});
  std::cout<<"Fid centered, now upsampling..."<<std::endl;
  cdouble r3 = cdouble(0.0, 1.0);
  //  wf_cen= dsp::hilbert_rot(wf_re, r,100);
  // std::vector<double> wf_po2 = dsp::hilbert(wf_pi_);
  cout<< "Fid centered, now upsampling"<<endl;
  // Make armadillo analytic version of wf_re 
  int M = 2*N;
  std::vector<double> wf_im = dsp::hilbert(wf_cen);*/
 
  int M = 2*N;
  arma::cx_vec wf(M);
  wf_cen = wf_fil;
  vector<double> wf_im = wf_fil_rot;
    
  // Artificially upsample.
  std::vector<double> wf_cen_up(M,0.0);
  std::vector<double> wf_im_up(M,0.0);
  auto it1 = wf_cen_up.begin();
  auto it2 = wf_im_up.begin();
  for (auto it3 = wf_cen.begin(); it3 !=wf_cen.end(); ++it3) {
    *(it1++) = *it3;
    *(it1++) = *it3;
  }

  for (auto it4 = wf_im.begin(); it4 !=wf_im.end(); ++it4) {
    *(it2++) = *it4;
    *(it2++) = *it4;
  }

  // Compare to taking the hilbert transform after the upsampling.
  cdouble r2 = cdouble(0.0 , 1.0);
  auto wf_im_test = dsp::hilbert_rot(wf_cen_up,r2);

  cout<<"Upsampling complete, now making analytic fid..."<<M<<endl;
  for (int i = 0; i < M; i++) {
    wf[i] = arma::cx_double(wf_cen_up[i], wf_im_up[i]);
  } 
  cout<<"Analytic fid ready, Now computing autocorrelation function"<<endl;

  // Window out the relevant portion of wf for speed efficiency.
  arma::cx_vec wf_temp(4*window);
  for (int i = 0; i< wf_temp.size(); i++) {
    wf_temp[i] = wf[1757-2*window +i];
  }
  cout<<" size of wf_temp is " << wf_temp.size()<<endl;
  
  // Now pad temp vector with zeros to provide more points for fitting.
  // wf_temp.resize(1500);

  //Compute auto-correlation function on relevant portion of data
  arma::cx_vec wf_rc = dsp::acorrelation(wf_temp,2*window  ,2*window);

  // Center real and imaginary acfs about zero.
  auto real = arma::conv_to<std::vector<double>>::from(arma::real(wf_rc));
  auto imag = arma::conv_to<std::vector<double>>::from(arma::imag(wf_rc));

  arma::vec Racf = arma::real(wf_rc);
  arma::vec Iacf = arma::imag(wf_rc);
  double rel_offset = (arma::min(Racf)-arma::min(Iacf));

  //  rel_offset = 0.0;
  
  Racf.transform([rel_offset](double val) {return val - rel_offset;});

  std::transform(real.begin(), real.end(),real.begin(),
                 [rel_offset](double val) {return val - rel_offset;});

  /* for (int i = 0; i < wf_rc.size(); i++) {
    wf_rc[i] = arma::cx_double(real[i], imag[i]);
    } */
  cout<< "relative offset is "<< rel_offset<< endl;
  wf_rc = arma::cx_vec(Racf, Iacf);
  auto wf_rc_beg = wf_rc.subvec(0,2*window -1);
  auto wf_rc_end = wf_rc.subvec(2*window,wf_rc.size()-1);
  arma::cx_vec offset (2000-2*window, arma::fill::zeros);
  arma::cx_vec union1 = arma::join_cols(wf_rc_beg, offset);
  wf_rc = arma::join_cols(union1, wf_rc_end);
 
  real = arma::conv_to<std::vector<double>>::from(arma::real(wf_rc));
  imag = arma::conv_to<std::vector<double>>::from(arma::imag(wf_rc));

  // Compute the power spectrum of the waveform.
  auto fft_wf = dsp::psd(wf_fil);
  auto fft_wf_po2 = dsp::psd(wf_im);
  auto acorr = arma::conv_to<std::vector<cdouble>>::from(wf_rc);
  auto fft_acorr = dsp::fft(acorr);//returns vector of size N/2
  std::vector<double> freqs_orig = dsp::fftfreq(wf_fil.size(), dt);
  std::vector<double> freqs = dsp::fftfreq(wf_rc.size(), dt);

  // Define vectors to draw on first part of canvas.

  TMultiGraph mg;
  TGraph gr = TGraph(3*window, &tm[fid_begin], &real[0]);gr.SetLineColor(kBlue);
  TGraph gi = TGraph(3*window, &tm[fid_begin], &imag[imag.size()-3*window]);gi.SetLineColor(kRed);
  TGraph psd = TGraph(freqs_orig.size()/2, &freqs_orig[1], &fft_wf[1]);psd.SetLineColor(kBlue);
  TGraph psd_quad = TGraph(freqs_orig.size()/2 , &freqs_orig[1], &fft_wf_po2[1]);psd_quad.SetLineColor(kRed);
  TCanvas c1; 
  c1.Divide(2,2);
  c1.cd(1);  
  mg.Add(&psd, "cp");
  mg.Add(&psd_quad, "cp");
  // mg.Add(&gi, "cp");
  //  mg.Add(&gi, "cp"); 

  TPad *p1 = (TPad *)(c1.cd(1));
  p1->SetLogy();

  mg.Draw("ap");

  cout<< "made it past the first part of the canvas"<< endl;

  //Define vectors to draw on the second part of the canvas.
  c1.cd(2);
  auto fft_acf = arma::fft(wf_rc);
  std::vector<double> real_fft = arma::conv_to<std::vector<double>>::from(arma::real(fft_acf));
  std::vector<double> imag_fft = arma::conv_to<std::vector<double>>::from(arma::imag(fft_acf));

  TMultiGraph mg1;
  TGraph gr1 = TGraph(wf_temp.size()/2, &freqs[0], &real_fft[0]);gr1.SetLineColor(kBlue);
  TGraph gi1 = TGraph(wf_temp.size()/2, &freqs[0], &imag_fft[0]);gi1.SetLineColor(kRed);

  mg1.Add(&gr, "cp");
  mg1.Add(&gi, "cp");

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
  std::vector<double> wf_test = dsp::hilbert_rot(wf_re, cdouble(0.0, 0.0),5);
  std::vector<double> wf_test_hil = dsp::hilbert_rot(wf_test, cdouble(0.0, 1.0));
  TMultiGraph mg2;
  TGraph gr2 = TGraph(wf_test.size()/5, &tm[fid_begin], &wf_test[0]);gr2.SetLineColor(kBlue);
  TGraph gi2 = TGraph(wf_test.size()/5, &tm[fid_begin], &wf_test_hil[0]);gi2.SetLineColor(kRed);

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

  for (int i =1 ; i<n; i++) {
    if (real_fft[max_bin+i]>.3*real_fft[max_bin]) {
      i_range = max_bin -1 - i;
      f_range = max_bin +1 + i;
      n = f_range-i_range;
    
    // if ((real_fft[max_bin-i]>.25*real_fft[max_bin])&(max_bin-i>0)) {
    //  i_range = max_bin - i;
      //   n=f_range-i_range;
    }else if (i_range <0){
      i_range = max_bin-2;
      f_range = max_bin+2;
    }
    //  n=f_range-i_range;
  }

  //Remove offset for better fit.
  int range = f_range - i_range;
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
  TF1 fit_func = TF1("fit_func", gaussian.c_str(), freqs[i_range],freqs[f_range]);
  fit_func.SetParameter(0, freqs[max_bin]);
  fit_func.SetParameter(1, 2);
  fit_func.SetParameter(2, peak.max());
  fit_func.SetParameter(3, 0);

  gr3.Fit(&fit_func, "R");

  std::cout<<"Center Frequency of FFT is: " << fit_func.GetParameter(0)<<endl;
  std::cout<<"Center Frequency from Centroid is: "<<my_fid.CalcCentroidFreq();

  //Save canvas.
  c1.Print("test_data/acorr_test_hil_up.pdf");

  return 0;
}
