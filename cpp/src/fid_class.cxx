#include "fid_class.h"

namespace fid {

FID::FID(const vec& wf, const vec& tm)
{
  // Copy the waveform and time to member vectors
  wf_ = wf;
  tm_ = tm;

  // Resize the temp array (maybe others?)
  temp_.reserve(wf_.size());

  // Initialize the FID for analysis
  CenterFid();
  CalcNoise();
  CalcMaxAmp();
  FindFidRange();
  CalcPowerEnvAndPhase();
  CalcFftFreq();
  GuessFitParams();
}

void FID::PrintDiagnosticInfo()
{
  cout << endl;
  cout << "Printing Diagostic Information for FID @ " << this << endl;
  cout << "noise level: " << noise_ << endl;
  cout << "waveform start, stop: " << i_wf_ << ", " << f_wf_ << endl;
  cout << "spectral start, stop: " << i_fft_ << ", " << f_fft_ << endl;
}

void FID::PrintDiagnosticInfo(std::iostream out)
{
  out << endl;
  out << "Printing Diagostic Information for FID @ " << this << endl;
  out << "noise level: " << noise_ << endl;
  out << "waveform start, stop: " << i_wf_ << ", " << f_wf_ << endl;
  out << "spectral start, stop: " << i_fft_ << ", " << f_fft_ << endl;
}

void FID::CenterFid()
{
  int w = params::zc_width;
  double sum  = std::accumulate(wf_.begin(), wf_.begin() + w, 0.0);
  double avg = sum / w; // to pass into lambda
  mean_ = avg; // save to class

  std::for_each(wf_.begin(), wf_.end(), [avg](double& x){ x -= avg; });
}

void FID::CalcNoise()
{ 
  // Grab a new handle to the noise window width for aesthetics.
  int w = params::zc_width;

  // Find the noise level in the head and tail.
  double head = stdev(wf_.begin(), wf_.begin() + w);
  double tail = stdev(wf_.rbegin(), wf_.rbegin() + w);

  // Take the smaller of the two.
  noise_ = (tail < head) ? (tail) : (head);
}

void FID::CalcMaxAmp() 
{
  auto mm = std::minmax_element(wf_.begin(), wf_.end());
  if (abs(*mm.first) > abs(*mm.second)) {
    max_amp_ = abs(*mm.first);
  } else {
    max_amp_ = abs(*mm.second);
  }
}

void FID::FindFidRange()
{
  // Find the starting and ending points
  double thresh = params::start_thresh * noise_;
  bool checks_out = false;

  // Find the first element with magnitude larger than thresh
  auto it_1 = wf_.begin();
  while (!checks_out) {

    auto it_i = std::find_if(it_1, wf_.end(), 
        [thresh](double x){return std::abs(x) > thresh;});

    if (it_i != wf_.end() && it_i+1 != wf_.end()) {
      checks_out = std::abs(*(it_i+1)) > thresh;
      it_1 = it_i + 1;

      // Turn the iterator into an index
      if (checks_out) {
        i_wf_ = std::distance(wf_.begin(), it_i);
      }

    } else {
        i_wf_ = std::distance(wf_.begin(), wf_.end());
      break;
    }
  }

  // Find the last element with magnitude larger than thresh
  checks_out = false;
  auto it_2 = wf_.rbegin();
  while (!checks_out) {

    auto it_f = std::find_if(it_2, wf_.rend(), 
      [thresh](double x){return std::abs(x) > thresh;});

    if (it_f != wf_.rend() && it_f+1 != wf_.rend()) {
      checks_out = std::abs(*(it_f+1)) > thresh;
      it_2 = it_f + 1;

      // Turn the iterator into an index
      if (checks_out) {
        f_wf_ = std::distance(it_f, wf_.rend());
      }

    } else {

      f_wf_ = std::distance(wf_.rend(), wf_.rend());
      break;
    }
  }

  // Mark the signal as bad if it didn't find signal above threshold.
  if (i_wf_ > wf_.size() * 0.9 || i_wf_ >= f_wf_) {

    isgood_ = false;
    i_wf_ = 0;
    f_wf_ = wf_.size() * 0.01;

  } else {

    isgood_ = true;
  }
}

void FID::CalcPowerEnvAndPhase()
{
  // Get the fft of the waveform first.
  auto fid_fft = dsp::fft(wf_);

  // Now get the imaginary harmonic complement to the waveform.
  auto wf_im = dsp::hilbert(fid_fft);

  // Optional lowpass filter. 
//  double df = (wf.size() - 1) / (wf.size() * (tm_[wf.size() - 1] - tm_[0]));
//  double cutoff_index = params::low_pass_freq / df;
//  

  // Now we can get power, envelope and phase.
  power_ = dsp::psd(fid_fft);
  phase_ = dsp::phase(wf_, wf_im);
  env_ = dsp::envelope(wf_, wf_im);
}

void FID::CalcFftFreq()
{
  // @todo: consider storing as start, step, stop
  freq_ = dsp::fftfreq(wf_);
}

void FID::GuessFitParams()
{
  // Guess the general fit parameters
  guess_.assign(5, 0.0);

  double f_mean;
  double f_mean2;
  double den;

  // find max index and set the fit window
  int max_idx = std::distance(power_.begin(),
    std::max_element(power_.begin(), power_.end()));

  if (max_idx - params::fit_width < 0) {
    i_fft_ = 0;

  } else {
    i_fft_ = max_idx - params::fit_width;
  }

  f_fft_ = max_idx + params::fit_width;
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
  freq_err_ = f_fit_.GetParError(0) / kTau;

  res_.resize(0);
  for (uint i = i_fft_; i < f_fft_ + 1; ++i){
    res_.push_back(power_[i] - func.Eval(freq_[i]));
  }

  return;
}

double FID::CalcZeroCountFreq()
{
  // set up vectors to hold relevant stuff about the important part
  temp_.resize(f_wf_ - i_wf_);

  int nzeros = 0;
  bool pos = wf_[i_wf_] >= 0.0;
  bool hyst = false;
  
  auto mm = std::minmax_element(wf_.begin(), wf_.end()); // returns {&min, &max}

  double max = *mm.second;
  if (std::abs(*mm.first) > max) max = std::abs(*mm.first);
  
  //  double max = (-(*mm.first) > *mm.second) ? -(*mm.first) : *mm.second;
  double thresh = params::hyst_thresh * max;

  int i_zero = -1;
  int f_zero = -1;

  // iterate over vector
  for (uint i = i_wf_; i < (f_wf_ - i_wf_); i++){

    // hysteresis check
    if (hyst){
      hyst = std::abs(wf_[i]) > thresh;
      continue;
    }

    // check for a sign change
    if ((wf_[i] >= 0.0) != pos){
      nzeros++;
      f_zero = i;
      if (i_zero == -1) i_zero = i;
      pos = !pos;
      hyst = true;
    }
  }

  // Use linear interpolation to find the zeros more accurately
  int i = i_zero;
  int f = f_zero;

  // do the interpolation
  double frac = std::abs(wf_[i] / (wf_[i-1] - wf_[i]));
  double ti = frac * tm_[i-1] + (1.0 - frac) * tm_[i];

  frac = std::abs(wf_[f] / (wf_[f-1] - wf_[f]));
  double tf = frac * tm_[f-1] + (1.0 - frac) * tm_[f];
  double freq = 0.5 * (nzeros - 1) / (tf - ti);

  // Calculate the error as done in the Priegl paper.
  freq_err_ = freq / ((tm_[f_wf_] - tm_[i_wf_]) * (noise_ * params::start_thresh));

  return freq;
}

double FID::CalcCentroidFreq()
{
  // Find the peak power
  double thresh = *std::max_element(power_.begin(), power_.end());
  thresh *= params::centroid_thresh;

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

  // Use the width of the peak, the stdev as the error estimator.
  freq_err_ = fid::stdev()
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
  for (uint i = 0; i < 5; i++){
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


double FID::CalcPhaseFreq(int poln)
{
  gr_time_series_ = TGraph(f_wf_ - i_wf_, &tm_[i_wf_], &phase_[i_wf_]);

  // Now set up the polynomial phase fit
  char fcn[20];
  sprintf(fcn, "pol%d", poln);
  f_fit_ = TF1("f_fit_", fcn);

  // Set the parameter guesses
  f_fit_.SetParameter(1, guess_[0] * kTau);

  // Adjust to ignore the edges
  int i = i_wf_ + params::edge_ignore;
  int f = f_wf_ - params::edge_ignore;

  // Do the fit.
  gr_time_series_.Fit(&f_fit_, "QMRSEX0", "C", tm_[i], tm_[f]);

  res_.resize(0);
  for (uint i = i_fft_; i < f_fft_ + 1; ++i){
    res_.push_back(wf_[i] - f_fit_.Eval(tm_[i]));
  }

  chi2_ = f_fit_.GetChisquare();
  freq_err_ = f_fit_.GetParError(1) / kTau;
  return f_fit_.GetParameter(1) / kTau;
}

double FID::CalcPhaseDerivFreq(int poln)
{
  // Do the normal polynomial phase fit to set the fit function
  CalcPhaseFreq(poln);

  // Find the initial phase by looking at the function's derivative
  return f_fit_.Derivative(tm_[i_wf_]) / kTau;
} 


double FID::CalcSinusoidFreq()
{
  // Normalize the waveform by the envelope
  temp_.resize(wf_.size());

  std::transform(wf_.begin(), wf_.end(), env_.begin(), temp_.begin(),
    [](double x_wf, double x_env) {return x_wf / x_env;});

  gr_time_series_ = TGraph(f_wf_ - i_wf_, &tm_[i_wf_], &temp_[i_wf_]);    

  f_fit_ = TF1("f_fit_", "[1] * sin([0] * x) + [2] * cos([0] * x)");

  // Guess parameters
  f_fit_.SetParameter(0, kTau * CalcZeroCountFreq());
  f_fit_.SetParameter(1, 1.0);
  f_fit_.SetParameter(2, phase_[i_wf_]);

  f_fit_.SetParLimits(1, 0.9, 1.1); // Should be exactly 1.0
  f_fit_.SetParLimits(2, -0.5 * kTau, 0.5 * kTau); // it's a phase

  // Adjust to ignore the edges
  int i = i_wf_ + params::edge_ignore;
  int f = f_wf_ - params::edge_ignore;

  // Do the fit.
  gr_time_series_.Fit(&f_fit_, "QMRSEX0", "C", tm_[i], tm_[f]);

  res_.resize(0);
  for (uint i = i_fft_; i < f_fft_ + 1; ++i){
    res_.push_back(wf_[i] - f_fit_.Eval(tm_[i]));
  }

  chi2_ = f_fit_.GetChisquare();
  freq_err_ = f_fit_.GetParError(0) / kTau;
  return f_fit_.GetParameter(0) / kTau;
}

} // fid
