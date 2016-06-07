#include "fid/wvd_fid.h"


namespace fid {

WvdFid::WvdFid(const std::string& fid_file, const int fid_stop, const int window, const bool upsample)
{
  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_file);
  
  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 
  fid_stop_ = fid_stop;

  WvdInit();
}


WvdFid::WvdFid(const char* fid_file, const int fid_stop, const int window, const bool upsample)
{
  // Convert the char pointer toa string.
  std::string fid_string(fid_file);

  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_string);
  
  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 
  fid_stop_ = fid_stop;

  WvdInit();
}


WvdFid::WvdFid(const std::vector<double>& wf, const std::vector<double>& tm, const int fid_stop, const int window, const bool upsample)
{
  // Copy the waveform and time to member vectors.
  wf_ = wf;
  tm_ = tm;

  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 
  fid_stop_ = fid_stop;

  WvdInit();
}

WvdFid::WvdFid(const std::vector<double>& wf, const int fid_stop, const int window, const bool upsample)
{
  // Copy the waveform and construct a generic time range.
  wf_ = wf;
  tm_ = construct_range(0.0, (double)wf_.size(), 1.0);
  
  // Set internal WVD parameters.
  upsample_ = upsample;
  window_ = window; 
  fid_stop_ = fid_stop;

  WvdInit();
}


void WvdFid::WvdInit()
{
  // Initialize relevant parameters.
  dt_ = std::abs(tm_[0]-tm_[1]);
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
  FindFidRange();
  AnalyticFid();

 // Calculate a health based on signal to noise.
  if (max_amp_ < noise_ * snr_thresh_) {
    health_ *= max_amp_ / (noise_ * snr_thresh_);
  }

  // And factor fid duration into health.
  if (f_wf_ - i_wf_ < wf_.size() * len_thresh_) {
    health_ *= (f_wf_ - i_wf_) / (wf_.size() * len_thresh_);
  }
}

void WvdFid::WvdCenter()
{
  // vector to store centered fid
  wf_cen_.resize(wf_.size(), 0.0);

  // hilbert transformation to produce centered fid.
  // cdouble r = cdouble(0.0 , 0.0);
  /*  std::vector<double> wf_pi = HilbertRot(wf_temp , r);
  std::transform( wf_.begin(), wf_.end(), wf_pi.begin(), wf_cen_.begin(),
                  [](double wf, double pi) {return wf - std::abs(wf - pi);});
  */      
  
  // The new way of creating a centered fid.
  // This Fid is shifted pi from the original and has a high pass filter applied so that 
  // successive hilbert transforms are better able to replicate the envelope.
  cdouble r = cdouble(0.0 , 0.0);
  wf_cen_ = HilbertRot(wf_, r, 10);
                                        
}

void  WvdFid::AnalyticFid()
{
  // Extract desired sinusoidal "FID" portion of the signal.
  unsigned int fid_begin = i_wf_;
  wf_cen_.erase(wf_cen_.begin(), wf_cen_.begin() + fid_begin);
  if ((fid_stop_ != 0)&(fid_stop_ +2*window_< wf_cen_.size())) wf_cen_.erase(wf_cen_.begin()+ fid_stop_+2*window_, wf_cen_.end());
    
  
  // Artificially double sampling rate
  int M, N;
  if (upsample_) {

    M = 2 * wf_cen_.size();
    N = wf_cen_.size();

  } else {

    M = wf_cen_.size();
    N = wf_cen_.size();
  }

  auto wf_im = HilbertRot(wf_cen_, cdouble(0.0,1.0));

  wf_cen_up_.resize(M, 0.0);
  std::vector<double> wf_im_up(M, 0.0);

  auto it1 = wf_cen_up_.begin();
  auto it2 = wf_im_up.begin();
  for (auto it3 = wf_cen_.begin(); it3 != wf_cen_.end(); ++it3) {
    *(it1++) = *it3;
    if (upsample_) {
      *(it1++) = *it3;
    }
  }
  
 for (auto it4 = wf_im.begin(); it4 != wf_im.end(); ++it4) {
    *(it2++) = *it4;
    if (upsample_) {
      *(it2++) = *it4;
    }
  }
  
  //Make signal analytic and centered about zero.
  wf_analytic_.resize(M, cdouble(0.0,0.0));

  for (uint i = 0; i < M; ++i) {
    wf_analytic_[i] = cdouble(wf_cen_up_[i], wf_im_up[i]);
  }
}

std::vector<cdouble> WvdFid::Filter(const std::vector<cdouble> &v, const int freq_cut)
{
  double G=0.0;
  double tmp=0.0;
  std::vector<cdouble> filter(v.size(), cdouble(0.0,0.0));
  int j = v.size();
  for (int i = 1; i<v.size(); i++) {
    // 4th Order Bessel High Pass Filter
     G = pow(freq_cut/(i/(v.size()*dt_)), 4)+ 
          10*pow(freq_cut/(2*i/(v.size()*dt_)), 3)+
          45*pow(freq_cut/(2*i/(v.size()*dt_)), 2)+
          105*pow(freq_cut/(2*i/(v.size()*dt_)), 1)+
          105;
     
     tmp =pow(1.0/G,0.5);
     filter[i] = cdouble(tmp,tmp);
  }
  
  return filter;
}

std::vector<double> WvdFid::HilbertRot(const std::vector<double>&v ,const cdouble i, const int freq_cut )
{
  //Return the call to the fft version.
  auto fft_vec = dsp::rfft(v);

  // Apply filter if called
  if (freq_cut !=0) {
    std::vector<cdouble> filter = Filter(fft_vec,freq_cut);
    std::transform(fft_vec.begin(), fft_vec.end(), filter.begin(), fft_vec.begin(),
                   [] (cdouble fft, cdouble fil) {return fft*fil;});
  }
  
  //Multiply by i^2.
  if (i != cdouble(0.0 , 0.0) ) {
    for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
      *it = i*(*it);
    }
  }
  // Reverse the fft.
  return dsp::irfft(fft_vec);
}

  std::vector<cdouble> WvdFid::AutoCorr( const std::vector<cdouble>& wf_temp, int idx)
{
  // Set Correlation Range
  int taumax = std::min({idx, window_});

  //make sure taumax is odd
  if ((taumax % 2 ==0) &(taumax !=0)) taumax -= 1;
  //  std::cout<< "Correlation window range is : " << taumax<<std::endl;

  // Set size of ACF vector
  std::vector<cdouble> acf(2*(taumax+1), cdouble(0.0,0.0));
  int N = 2*(taumax+1);

  // Construct auto-correlation function that gives purely real fft.
  // Also optimized to produce symmetric fft peak.
  for (int i =-taumax+1; i < taumax ; i++) {
    if ((i < 0)&(idx+taumax<N)) {
      acf[N+i] = std::conj(wf_temp[idx+i])*wf_temp[idx-i];//starts filling acf at index 0
    }   
    if ((i >= 0)&(idx+taumax<N)) {
      acf[i] = std::conj(wf_temp[idx+i])*wf_temp[idx-i];
    }
  }
 
  return acf;
} 

double WvdFid::WvdZeroFreq(const int fid_stop)
{
  // Make sure FID is healthy.
  if (health_>99) {
    // Extract Frequency from initial part of the fid.
    fid_stop_ = fid_stop;
  
    // Now Extract frequencies.
    WvdFreqExt();

    // Make a Histogram of extracted frequencies.
    auto extrema = std::minmax_element(wvd_f_.begin()+100, wvd_f_.end()-50);
    TH1D h1 = TH1D("h1", "Extracted Frequency Distribution", 
                   1000,*extrema.first, *extrema.second);
    for (int i = 0; i< wvd_f_.size(); i++) {
      h1.Fill(wvd_f_[i]);
    }

    // Get frequency corresponding to max bin.
    int maxbin = h1.GetMaximumBin();
    wvd_f0_ = h1.GetXaxis()->GetBinCenter(maxbin);
  } else {
    wvd_f0_ = 0.0;
  }

  return wvd_f0_; 
}

void WvdFid::WvdFreqExt() 
{
  if (health_>99) {
    // Set size of the vector holding extracted frequencies.
    int N = wf_cen_.size();
    if (fid_stop_ ==0) {
      wvd_f_.resize(N, 0.0);
    }else{
      wvd_f_.resize(fid_stop_, 0.0);
    }

    // @Todo:guess optimal WVD parameters
    int n_temp = 2000;

    // Get frequency dimentionality 
    std::vector<double> freq = dsp::fftfreq(n_temp, dt_);
  
    // Now extract frequencies for each time instant using psuedo Wigner-Ville distribution 
    // I need to determine a better way of determing the start point.
    // 50 was chosen from quick crash tests.  
    for (auto it = wvd_f_.begin()+50; it != wvd_f_.end()-50; ++it) {

      // Define some relevant quantities
      int idx = std::distance(wvd_f_.begin(), it);
      std::vector<cdouble> wf_temp(2*window_, cdouble(0.0,0.0));
      arma::cx_vec auto_corr(wf_temp.size());
      auto_corr.zeros();

      //  std::cout<<"FFT for loop iterator has a value of: "<<idx<<std::endl;

      // Create smaller windowed temp version of wf_analytic_ for efficiency
      // and compute the autocorrelation function on smaller temp waveform.
      if (idx < window_) {
        for (int i = 0; i< wf_temp.size(); i++) {
          wf_temp[i] = wf_analytic_[i];
        }
 
        auto_corr= arma::conv_to<arma::cx_vec>::from(AutoCorr(wf_temp, idx));
      
      } else if (idx > (N-1-window_)) {
        for (int i = 0; i< wf_temp.size(); i++) {
          wf_temp[i] = wf_analytic_[N-1-wf_temp.size()+i];
        }
  
        auto_corr= arma::conv_to<arma::cx_vec>::from(AutoCorr(wf_temp,N-1-idx));
      
      } else { 
        for (int i = 0; i< wf_temp.size(); i++) {
          wf_temp[i] = wf_analytic_[idx-window_+i+1];
        }
   
        auto_corr= arma::conv_to<arma::cx_vec>::from(AutoCorr(wf_temp, window_));
       
      }
    
      // Adjust offset between quadrature components of acf. 
      arma::vec Racf = arma::real(auto_corr);
      arma::vec Iacf = arma::imag(auto_corr);
      double rel_offset = (arma::min(Racf)-arma::min(Iacf))/2;
    
      Racf.transform([rel_offset](double val) {return val - rel_offset;});
      auto_corr = arma::cx_vec(Racf, Iacf);

      // Pad auto_corr with zeros in the middle for better fit results
      int n_acf = auto_corr.size();
      auto acbegin = auto_corr.subvec(0,n_acf/2-1);
      auto acend = auto_corr.subvec(n_acf/2,n_acf-1);
      arma::cx_vec zeros_ (n_temp - n_acf, arma::fill::zeros);

      arma::cx_vec acbegin_zeros = arma::join_cols(acbegin, zeros_);
      auto_corr = arma::join_cols(acbegin_zeros, acend);
      acf_ = arma::conv_to<std::vector<cdouble>>::from(auto_corr);
   
      // Compute the FFT of the AutoCorrelation function
      arma::vec re_fft_acorr = arma::real(arma::fft(auto_corr));

      // Find Peak bin and range about peak to extract FFT frequency.
      arma::uword max_bin;
      double max_val = re_fft_acorr.max(max_bin);
   
      int i_range=0; int f_range=0;
      int n = re_fft_acorr.size();
      for (int i =1 ; i<n; i++) {
        if (re_fft_acorr[max_bin+i]>.3*re_fft_acorr[max_bin]) {
          i_range = max_bin- 1 - i;
          f_range = max_bin+ 1 + i;
          n = f_range-i_range;
        } else if (i_range <0){
          i_range = max_bin-2;
          f_range = max_bin+2;
        }
      }
    
      //create vector of normalized peak for better fitting.
      int range = f_range - i_range;
      arma::vec peak = arma::zeros(range);
      double sum = 0;
      double dx = freq[2]-freq[1];
      for (int i=i_range; i<f_range+1; i++) {
        sum += .5*(re_fft_acorr[i+1]+re_fft_acorr[i])*dx;
      }
  
      for (int i = 0; i<peak.size(); i++) {
        peak[i]=re_fft_acorr[i_range+i]/sum;
      }
    
      // Now perform fit to Gaussian
      std::string gaussian("[2]*exp(-(x-[0])^2/(2*[1]^2))+[3]");
      TGraph gr3 = TGraph(range, &freq[i_range], &peak[0]);
      TF1 fit_func = TF1("fit_func", gaussian.c_str(), freq[i_range],freq[f_range]);
  
      fit_func.SetParameter(0, freq[max_bin]);
      fit_func.SetParameter(1, 2);
      fit_func.SetParameter(2, peak.max());
      fit_func.SetParameter(3, 0);
    
      gr3.Fit(&fit_func, "RQ");

      // Save the extracted frequency
      wvd_f_[idx] = fit_func.GetParameter(0);

    } 
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
