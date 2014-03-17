#include "fid_class.h"

namespace fid{

  FID::FID(const vec& wf, const vec& tm)
  {
    // Copy the waveform and time to member vectors
    wf_ = wf;
    tm_ = tm;

    // Resize the temp array (maybe others?)
    temp_.reserve(wf_.size());

    // Initialize the FID for analysis
    CalcNoise();
    FindFidRange();
    CalcPowerEnvAndPhase();
    CalcFftFreq();
    GuessFitParams();
  }

  FID::~FID(){
    // TODO
  }

    
  void FID::CalcNoise(){
    // Find the noise level

    double noise = std::accumulate(wf_.begin(), wf_.begin() + kZCWidth, 0.0, 
                    [](double x, double y) {return x + y * y;}); 
    double temp = std::accumulate(wf_.rbegin(), wf_.rbegin() + kZCWidth, 0.0,
                    [](double x, double y) {return x + y * y;});

    noise = (temp < noise) ? (temp / kZCWidth) : (noise / kZCWidth);

    // Set the member noise
    noise_ = std::pow(noise, 0.5); // Want the RMS
  }


  void FID::FindFidRange(){
    // Find the starting and ending points

    double thresh = kStartThresh * noise_;

    // Find the first element with magnitude larger than thresh
    auto it_i = std::find_if(wf_.begin(), wf_.end(), 
        [thresh](double x){return std::abs(x) > thresh;});

    // Find the last element with magnitude larger than thresh
    auto it_f = std::find_if(wf_.rbegin(), wf_.rend(), 
        [thresh](double x){return std::abs(x) > thresh;});

    // Turn the iterators into indexes
    i_wf_ = std::distance(wf_.begin(), it_i);
    f_wf_ = std::distance(it_f, wf_.rend());
  }


  void FID::CalcPowerEnvAndPhase(){
    // Do the FFT and find the power and envelope function

    // Set up the FFT plan
    int N = wf_.size();  // size of data
    int n = N / 2 + 1;      // size of rfft
    double Nroot = std::sqrt(N);
    double temp = 0.0;

    // Resize vectors
    freq_.resize(n);
    power_.resize(n);
    phase_.resize(N);

    // To store inverse transform
    vec wf_im(N, 0.0);

    // Plan and execute the FFT
    fftw_complex *fft = new fftw_complex[n];
    fftw_complex *ffti = new fftw_complex[n];

    fftw_plan wf_to_fft;
    wf_to_fft = fftw_plan_dft_r2c_1d(N, &wf_[0], &fft[0], FFTW_ESTIMATE);
    fftw_execute(wf_to_fft);

    // Store the power
    for (int i = 0; i < n; i++){
      power_[i] = std::pow(fft[i][0], 2) + std::pow(fft[i][1], 2);
    }

    // Apply low-pass Butterworth Filter
    double df = (N - 1) / (N * (tm_[N - 1] - tm_[0]));
    // convert cutoff frequency to fft index
    double cutoff_index = kLowPassFreq / df;

    // Apply low pass and get fft * (-i)
    for (int i = 0; i < n; i++){
      fft[i][0] *= LowPassFilter(i, cutoff_index, n) / Nroot;
      fft[i][1] *= LowPassFilter(i, cutoff_index, n) / Nroot;
      ffti[i][0] = fft[i][1];
      ffti[i][1] = -1 * fft[i][0];
    }

    // Now get the Hilbert transform
    fftw_plan fft_to_wf;
    fft_to_wf = fftw_plan_dft_c2r_1d(N, &fft[0], &wf_im[0], FFTW_ESTIMATE);
    fftw_execute(fft_to_wf);

    // fftw is unnormalized fft, so rescale
    for (auto it = wf_im.begin(); it != wf_im.end(); it++){
      *it /= Nroot;
    }

    // Set the envelope function
    env_.resize(N);
    std::transform(wf_.begin(), wf_.end(), wf_im.begin(), env_.begin(),
      [](double re, double im) {return std::sqrt(im * im + re * re);});

    // Store the moduloed phase
    phase_.resize(N);
    std::transform(wf_.begin(), wf_.end(), wf_im.begin(), phase_.begin(),
      [](double re, double im) {return std::atan2(im, re);});

    // Now unwrap the phase
    double thresh = kMaxPhaseJump;
    for (auto it = phase_.begin() + i_wf_ + 1; it != phase_.end(); it++){
      static int k;
      if (it == phase_.begin() + i_wf_ + 1) k = 0; // reset

      // Add current total
      *it += k * kTau;

      // Check for jumps, both positive and negative.
      if (*(it-1) - *it > thresh){
        k++;
        *it += kTau;

      } else if (*it - *(it-1) > thresh){
        k--;
        *it -= kTau;

      }
    }

    // Clean up
    delete[] fft;
    fftw_destroy_plan(wf_to_fft);
    fftw_destroy_plan(fft_to_wf);
  }

  void FID::CalcFftFreq(){
    // Set up an array of FFT frequencies
    // @bug This could be stored as first, last, step_size...
    int N = wf_.size();
    double df = (N - 1) / (N * (tm_[N - 1] - tm_[0]));

    if (N % 2 == 0){

      freq_.reserve(N / 2 + 1);
      freq_.resize(0);
      
      for (int i = 0; i < N / 2 + 1; i++){
        freq_.push_back(i * df);
      }

    } else {

      freq_.reserve((N + 1) / 2);
      freq_.resize(0);

      for (int i = 0; i < (N + 1) / 2 + 1; i++){
        freq_.push_back(i * df);
      }
    }

    return;
  }

  void FID::GuessFitParams(){
    // Guess the general fit parameters
    guess_.assign(5, 0.0);

    double f_mean;
    double f_mean2;
    double den;

    // find max index and set the fit window
    int max_idx = std::distance(power_.begin(),
      std::max_element(power_.begin(), power_.end()));

    i_fft_ = max_idx - kFitWidth;
    if (i_fft_ < 0) i_fft_ = 0;

    f_fft_ = max_idx + kFitWidth;
    if (f_fft_ > power_.size()) f_fft_ = power_.size();

    auto it_pi = power_.begin() + i_fft_; // to shorten subsequent lines
    auto it_pf = power_.begin() + f_fft_;
    auto it_fi = freq_.begin() + i_fft_;

    // Compute some moments
    f_mean = std::inner_product(it_pi, it_pf, it_fi, 0.0);
    den = std::accumulate(it_pi, it_pf, 0.0); // normalization
    f_mean /= den;

    // find average power squared
    f_mean2 = std::inner_product(it_pi, it_pf, it_fi, 0.0,
      [](double sum, double x) {return sum + x;},
      [](double x1, double x2) {return x1 * x2 * x2;});
    f_mean2 /= den;

    // frequency
    guess_[0] = f_mean;

    // peak width
    guess_[1] = std::sqrt(f_mean2 - f_mean * f_mean);

    // amplitude
    guess_[2] = power_[max_idx];

    // background
    guess_[3] = noise_;

    // exponent
    guess_[4] = 2.0;

    return;
  }

  void FID::FreqFit(TF1& func)
  {
    // Make a TGraph to fit
    gr_freq_series_ = TGraph(f_fft_ - i_fft_, &freq_[i_fft_], &power_[i_fft_]);

    gr_freq_series_.Fit(&func, "QMRSEX0", "C", freq_[i_fft_], freq_[f_fft_]);

    chi2_ = func.GetChisquare();

    return;
  }

  double FID::CalcZeroCountFreq()
  {
    // set up vectors to hold relevant stuff about the important part
    temp_.resize(f_wf_ - i_wf_);
    vector<int> sign (f_wf_ - i_wf_, 0);
    vector<int> diff (f_wf_ - i_wf_, 0);

    // // smooth the waveform, exponential moving average algorithm
    // double a = kZCAlpha;
    // std::partial_sum(wf_.begin() + i_wf_, wf_.begin() + f_wf_, temp_.begin(), 
    //   [a](double sum, double x) {return sum * (1.0 - a) + x * a;});

    // get the sign
    std::transform(temp_.begin(), temp_.end(), sign.begin(),
      [](double x) {return x >= 0.0 ? 1 : -1;});

    // Find out when the sign changes
    std::adjacent_difference(sign.begin(), sign.end(), diff.begin());

    // Count the sign changes
    int zeros = std::count_if(diff.begin() + 1, diff.end(), 
      [](int x) {return x != 0;});

    // Use linear interpolation to find the zeros more accurately
    auto it_i = std::find_if(diff.begin() + 1, diff.end(),
      [](int x) {return x != 0;});

    auto it_f = std::find_if(diff.rbegin(), diff.rend(),
      [](int x) {return x != 0;});

    // convert iterators to indices
    int i = i_wf_ + std::distance(diff.begin(), it_i) - 1;
    int f = i_wf_ + std::distance(it_f, diff.rend()) - 1;

    // do the interpolation
    double frac = std::abs(wf_[i] / (wf_[i-1] - wf_[i]));
    double ti = frac * tm_[i-1] + (1.0 - frac) * tm_[i];

    frac = std::abs(wf_[f] / (wf_[f-1] - wf_[f]));
    double tf = frac * tm_[f-1] + (1.0 - frac) * tm_[f];

    return 0.5 * (zeros - 1) / (tf - ti);
  }

  double FID::CalcCentroidFreq()
  {
    // Find the peak power
    double thresh = *std::max_element(power_.begin(), power_.end());
    thresh *= kCentroidThresh;

    // Find the indices for a window around the max
    int it_i = std::distance(power_.begin(), 
      std::find_if(power_.begin(), power_.end(), 
        [thresh](double x) {return x > thresh;}));

    // reverse the iterators
    int it_f = -1 * std::distance(power_.rend(),
      std::find_if(power_.rbegin(), power_.rend(), 
        [thresh](double x) {return x > thresh;}));

    // Now compute the power weighted average
    double pwfreq = 0.0;
    double pwsum = 0.0;

    for (int i = it_i; i < it_f; i++){
      pwfreq += power_[i] * freq_[i];
      pwsum  += power_[i];
    }

    return pwfreq / pwsum;
  }


  double FID::CalcAnalyticalFreq()
  {
    // @todo check the algebra on this guy
    // Set the fit function
    std::string num1("[2] * ([0]^2 - 0.5 * [1] * [0] * sin(2 *[4])");
    std::string num2(" + (0.5 * [1])^2 + x^2 - [0]^2 * sin([4])^2) / ");
    std::string den1("((0.5 * [1])^2 + 2 * (x^2 - [0]^2) * (x^2 + [0]^2)");
    std::string den2(" + (x^2 - [0]^2)^2) + [3]");

    f_fit_ = TF1("f_fit_", (num1 + num2 + den1 + den2).c_str());

    // Set the parameter guesses
    for (int i = 0; i < 5; i++){
      f_fit_.SetParameter(i, guess_[i]);
    }

    // Limits
    f_fit_.SetParLimits(4, 0.0, kTau);

    FreqFit(f_fit_);
    return f_fit_.GetParameter(0);
  }

  
  double FID::CalcLorentzianFreq()
  {
    // Set the fit function
    f_fit_ = TF1("f_fit_", "[2] / (1 + ((x - [0]) / (0.5 * [1]))^2) + [3]");

    // Set the parameter guesses
    for (int i = 0; i < 4; i++){
      f_fit_.SetParameter(i, guess_[i]);
    }

    FreqFit(f_fit_);
    return f_fit_.GetParameter(0);
  }

  double FID::CalcSoftLorentzianFreq()
  {
    // Set the fit function
    f_fit_ = TF1("f_fit_", "[2] / (1 + ((x - [0]) / (0.5 * [1]))^[4]) + [3]");

    // Set the parameter guesses
    for (int i = 0; i < 5; i++){
      f_fit_.SetParameter(i, guess_[i]);
    }

    f_fit_.SetParLimits(4, 1.0, 3.0);

    FreqFit(f_fit_);
    return f_fit_.GetParameter(0);
  }

  double FID::CalcExponentialFreq()
  {
    // Set the fit function
    f_fit_ = TF1("f_fit_", "[2] * exp(-abs(x - [0]) / [1]) + [3]");

    // Set the parameter guesses
    for (int i = 0; i < 4; i++){
      f_fit_.SetParameter(i, guess_[i]);
    }

    f_fit_.SetParameter(1, 0.5 * guess_[1] / std::log(2));

    FreqFit(f_fit_);
    return f_fit_.GetParameter(0); 
  }

  double FID::CalcPhaseFreq(int n)
  {
    // 
    gr_time_series_ = TGraph(f_wf_ - i_wf_, &tm_[i_wf_], &phase_[i_wf_]);

    // Now set up the polynomial phase fit
    char fcn[20];
    sprintf(fcn, "pol%d", n);
    f_fit_ = TF1("f_fit_", fcn);

    // Set the parameter guesses
    f_fit_.SetParameter(1, guess_[0] * kTau);

    // Adjust to ignore the edges
    int i = i_wf_ + kEdgeIgnore;
    int f = f_wf_ - kEdgeIgnore;

    // Do the fit.
    gr_time_series_.Fit(&f_fit_, "QMRSEX0", "C", tm_[i], tm_[f]);

    chi2_ = f_fit_.GetChisquare();
    return f_fit_.GetParameter(1) / kTau;
  }

  double FID::CalcSinusoidFreq()
  {
    // Normalize the waveform by the envelope
    temp_.resize(wf_.size());

    std::transform(wf_.begin(), wf_.end(), env_.begin(), temp_.begin(),
      [](double x_wf, double x_env) {return x_wf / x_env;});

    gr_time_series_ = TGraph(f_wf_ - i_wf_, &tm_[i_wf_], &temp_[i_wf_]);    

    f_fit_ = TF1("f_fit_", "[1] * sin([0] * x + [2])");

    // Set guess parameters
    f_fit_.SetParameter(0, 23.0 * kTau);
    f_fit_.SetParameter(1, 1.0);
    f_fit_.SetParameter(2, phase_[i_wf_]);

    f_fit_.SetParLimits(1, 0.9, 1.1); // Should be exactly 1.0
    f_fit_.SetParLimits(2, -0.5 * kTau, 0.5 * kTau); // it's a phase

    // Adjust to ignore the edges
    int i = i_wf_ + kEdgeIgnore;
    int f = f_wf_ - kEdgeIgnore;

    // Do the fit.
    gr_time_series_.Fit(&f_fit_, "QMRSEX0", "C", tm_[i], tm_[f]);

    chi2_ = f_fit_.GetChisquare();
    return f_fit_.GetParameter(0) / kTau;
  }


} // fid
