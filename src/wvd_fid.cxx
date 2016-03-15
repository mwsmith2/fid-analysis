#include "fid/wvd_fid.h"

namespace fid {

WvdFid::WvdFid(const std::string& fid_file, const bool upsample, const int window)
{
  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_file);
  
  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 

  WvdInit();
}


WvdFid::WvdFid(const char* fid_file, const bool upsample, const int window)
{
  // Convert the char pointer toa string.
  std::string fid_string(fid_file);

  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_string);
  
  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 

  WvdInit();
}


WvdFid::WvdFid(const std::vector<double>& wf, const std::vector<double>& tm, const bool upsample, const int window)
{
  // Copy the waveform and time to member vectors.
  wf_ = wf;
  tm_ = tm;

  std::cout << "size of wf is "<< wf.size()<<std::endl;

  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 

  WvdInit();
}

WvdFid::WvdFid(const std::vector<double>& wf, const bool upsample, const int window)
{
  // Copy the waveform and construct a generic time range.
  wf_ = wf;
  tm_ = construct_range(0.0, (double)wf_.size(), 1.0);
  
  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 
  
  WvdInit();
} 

void WvdFid::WvdInit()
{
  // Initialize the health properly.
  // @Todo: properly implement health functionality
  health_ = 100.0;

  // Resize the temp array (maybe others?)
  temp_.reserve(wf_.size());
  wf_cen_.reserve(wf_.size());

  // Initialize the BaseFid for analysis
  LoadParams();
  CenterFid();
  CalcNoise();
  CalcMaxAmp();
  WvdCenter();
  AnalyticFid();

  // Default WVD frequency extraction technique.
  WvdFreqExt();
}

void WvdFid::WvdCenter()
{
  // vector to store centered fid
  wf_cen_.resize(wf_.size(), 0.0);

  // hilbert transformation to produce centered fid.
  cdouble r = cdouble(1.0 , 0.0);
  std::vector<double> wf_temp = wf_;
  std::vector<double> wf_pi = HilbertRot(wf_temp , r);
  std::transform( wf_.begin(), wf_.end(), wf_pi.begin(), wf_cen_.begin(),
                  [](double wf, double pi) {return wf - std::abs(wf - pi);});
  
}

void  WvdFid::AnalyticFid()
{
  // Extract Sinusoidal "FID" portion of the signal.
  FindFidRange();
  unsigned int fid_begin = i_wf_;
  wf_cen_.erase(wf_cen_.begin(), wf_cen_.begin() + fid_begin);
  
   //Artificially double sampling rate if desired
  int M, N;
  if (upsample_) {

    M = 2 * wf_cen_.size();
    N = wf_cen_.size();

  } else {

    M = wf_cen_.size();
    N = wf_cen_.size();
  }

  wf_cen_up_.resize(M, 0.0);

  auto it1 = wf_cen_up_.begin();
  for (auto it2 = wf_cen_.begin(); it2 != wf_cen_.end(); ++it2) {
    *(it1++) = *it2;
    if (upsample_) {
      *(it1++) = *it2;
    }
  }

  //Make signal analytic and centered about zero.
  wf_analytic_.resize(M, cdouble(0.0,0.0));

  auto wf_im = dsp::hilbert(wf_cen_up_);

  for (uint i = 0; i < M; ++i) {
    wf_analytic_[i] = cdouble(wf_cen_up_[i], wf_im[i]);
  }
}

std::vector<double> WvdFid::HilbertRot(const std::vector<double>&v ,const cdouble i)
{
  //Return the call to the fft version.
  auto fft_vec = dsp::rfft(v);

  //Multiply by i^2.
  for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
    if (it ==fft_vec.begin()) *it = cdouble(0,0);
    *it = i*(*it);
  }

  // Reverse the fft.
  return dsp::irfft(fft_vec);
}

std::vector<cdouble> WvdFid::AutoCorr(int idx)
{
  // Resize relevant vectors
  int N = wf_analytic_.size();
  acf_.resize(wf_analytic_.size(), cdouble(0.0,0.0));

  // Set Correlation Range
  int taumax = std::min({idx, std::abs(N/2-idx), window_});
  //make sure taumax is odd
  if ((taumax % 2 ==0) &(taumax !=0)) taumax -= 1;
  std::cout<< "Correlation window range is : " << taumax<<std::endl;

  // Construct auto-correlation function that gives purely real fft.
  // Also optimized to produce symmetric fft peak.
  for (int i =-taumax+1; i < taumax ; i++) {
    //    std::cout<< "correlation variable is " <<i <<" idx is "<< idx<<std::endl;
    if ((i < 0)&(idx+taumax<N)) {
      acf_[N+i] = std::conj(wf_analytic_[idx-i])*wf_analytic_[idx+i];//starts filling acf at index 0
    }   
    if ((i >= 0)&(idx+taumax<N)) {
      acf_[i] = std::conj(wf_analytic_[idx-i])*wf_analytic_[idx+i];
    }
  }
 
  //even odd symmetry test
  /* double even=0.0, odd=0.0, sum=0.0, prev =0.0;
  for (int i =0; i< N/2; i++) {
    even += (acf(i).real()-acf((N-1)-(i)).real());
    odd += (acf(i).imag()+acf((N-1)-(i)).imag());
    sum += abs(acf(i).imag()); 
  }
  */

  return acf_;
} 

void WvdFid::WvdFreqExt() 
{
  // Set size of vector
  int N = wf_cen_.size();
  wvd_f_.resize(N, 0.0);

  std::cout<< "size of wf_cen_ is " << wf_cen_.size()<<std::endl;
  
  // Now compute the Wigner-Ville Distribution
  // I need to determine a better way of determing the start point.
  // 10 was chosen from quick crash tests.

   // Get frequency dimenstionality 
  std::vector<double> freq = dsp::fftfreq(N, .001);
  
  for (int idx = i_wf_ + 10; idx < N-10; ++idx) {
    // Compute AutoCorrelation function for each time instant, idx.
    arma::cx_vec auto_corr= arma::conv_to<arma::cx_vec>::from(AutoCorr(idx));
    
    // Compute the FFT of the AutoCorrelation function
    arma::vec re_fft_acorr = arma::real(arma::fft(auto_corr));

    // Find Peak bin and range about peak to extract FFT frequency.
    arma::uword max_bin;
    double max_val = re_fft_acorr.max(max_bin);
    std::cout<< "peak bin via armadillo is: "<<max_bin<<std::endl;
   
    int i_range=0; int f_range = 0;
    int n = re_fft_acorr.size();
    for (int i =1 ; i<n+1; i++) {
      if (re_fft_acorr[max_bin+i]>.1*re_fft_acorr[max_bin]) {
        f_range= max_bin + i+1;
        n = f_range-i_range;
      }
      if ((max_bin-i>0)&(re_fft_acorr[max_bin-i]>.1*re_fft_acorr[max_bin])) {
        i_range = max_bin - i;
      }else if (i_range <0){
        i_range = 0;
        f_range = 2*max_bin;
      }
    }
    
    //create vector of normalized peak for better fitting.
    double range = f_range - i_range;
    arma::vec peak = arma::zeros(range);
    double sum = 0;
    double dx = freq[1]-freq[0];
    for (int i=i_range; i<f_range+1; i++) {
      sum += .5*(re_fft_acorr[i+1]+re_fft_acorr[i])*dx;
    }
  
    for (int i = 0; i<peak.size(); i++) {
      peak[i]=re_fft_acorr[i_range+i]/sum;
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

    std::cout<<"Center Frequency of FFT is: " << fit_func.GetParameter(0)/2<<std::endl;
    
    // Save the extracted frequency
    wvd_f_[idx] = fit_func.GetParameter(0)/2;

  }
}

void WvdFid::SaveWvd(std::string filename, std::string title)
{
   // If no title supplied give a reasonable default.
  if (title == "") {

    title = std::string("WVD; time [ms]; amplitude [a.u.]");

  } else {

    // In case they didn't append x/y labels.
    title.append("; time [ms]; amplitude [a.u.]");
  }
  TCanvas c1;
  c1.Print(filename.c_str());

  TGraph gr = TGraph(wvd_f_.size(), &tm_[i_wf_], &wvd_f_[0]);

  gr.SetTitle(title.c_str());
  gr.Draw("cp");
  
}

void WvdFid::MutltiPlot(const int pads)
{
  TCanvas c1; 
  for (int i = 0 ; i< pads ; i++) {
    
  }
}
} // fid
