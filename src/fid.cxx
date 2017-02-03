#include "fid/fid.h"
#include "fid/params_extdef.h"

namespace fid {

Fid::Fid(const std::string& fid_file)
{
  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_file);

  Init();
}

Fid::Fid(const char* fid_file)
{
  // Convert the char pointer toa string.
  std::string fid_string(fid_file);

  // Read and store the waveform and time from a .fid file.
  LoadTextData(fid_string);

  Init();
}


Fid::Fid(const std::vector<double>& wf, const std::vector<double>& tm)
{
  // Copy the waveform and time to member vectors.
  wf_ = wf;
  tm_ = tm;

  Init();
}

Fid::Fid(const std::vector<double>& wf)
{
  // Copy the waveform and construct a generic time range.
  wf_ = wf;
  tm_ = construct_range(0.0, (double)wf_.size(), 1.0);

  Init();
}


void Fid::InitHook()
{
  CalcPowerEnvAndPhase();
  CalcFftFreq();
  GuessFitParams();
}


double Fid::GetFreq(const std::string& method_name)
{
  return GetFreq(ParseMethod(method_name));
}


double Fid::GetFreq(const char* method_name)
{
  return GetFreq((std::string(method_name)));
}


double Fid::GetFreq(const Method m)
{
  // Recalculate if necessary
  if (freq_method_ != m) {

    freq_method_ = m;
    return CalcFreq();

  } else {

    return freq_;
  }
}


// Calculate the frequency using the current Method
double Fid::CalcFreq()
{
  switch(freq_method_) {

    case ZC:
    case ZEROCOUNT:
      CalcZeroCountFreq();
      break;

    case CN:
    case CENTROID:
      CalcCentroidFreq();
      break;

    case AN:
    case ANALYTICAL:
      CalcAnalyticalFreq();
      break;

    case LZ:
    case LORENTZIAN:
      CalcLorentzianFreq();
      break;

    case EX:
    case EXPONENTIAL:
      CalcExponentialFreq();
      break;

    case PH:
    case PHASE:
      CalcPhaseFreq();
      break;

    case SN:
    case SINUSOID:
      CalcSinusoidFreq();
      break;

    default:
      // Reset the Method to a valid one.
      freq_method_ = params::freq_method;
      CalcPhaseFreq();
      break;
  }

  return freq_;
}


void Fid::CalcPowerEnvAndPhase()
{
  // Get the fft of the waveform first.
  auto fid_fft = dsp::rfft(wf_);

  // Now get the imaginary harmonic complement to the waveform.
  auto wf_im = dsp::hilbert(wf_);

  // Optional lowpass filter. 
  //  double df = (wf.size() - 1) / (wf.size() * (tm_[wf.size() - 1] - tm_[0]));
  //  double cutoff_index = low_pass_freq / df;
  //  

  // Now we can get power, envelope and phase.
  psd_ = dsp::norm(fid_fft);
  phi_ = dsp::phase(wf_, wf_im);
  env_ = dsp::envelope(wf_, wf_im);
}

void Fid::CalcFftFreq()
{
  // @todo: consider storing as start, step, stop
  fftfreq_ = dsp::fftfreq(tm_);
}

void Fid::GuessFitParams()
{
  // Guess the general fit parameters
  guess_.assign(6, 0.0);

  double f_mean;
  double f_mean2;
  double den;

  // find max index and set the fit window
  int max_idx = std::distance(psd_.begin(),
    std::max_element(psd_.begin()+1, psd_.end()));

  i_fft_ = max_idx - fft_peak_width_;
  if (max_idx - fft_peak_width_ < 0) i_fft_ = 0;  

  f_fft_ = max_idx + fft_peak_width_;
  if (f_fft_ > psd_.size()) f_fft_ = psd_.size();

  auto it_pi = psd_.begin() + i_fft_; // to shorten subsequent lines
  auto it_pf = psd_.begin() + f_fft_;
  auto it_fi = fftfreq_.begin() + i_fft_;

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
  guess_[2] = psd_[max_idx];

  // background
  guess_[3] = noise_;

  // exponent
  guess_[4] = 2.0;

  return;
}


void Fid::FreqFit(TF1& func)
{
  // Make a TGraph to fit
  gr_freq_series_ = TGraph(f_fft_ - i_fft_, &fftfreq_[i_fft_], &psd_[i_fft_]);

  gr_freq_series_.Fit(&func, "QMRSEX0", "C", fftfreq_[i_fft_], fftfreq_[f_fft_]);

  chi2_ = func.GetChisquare();
  freq_err_ = f_fit_.GetParError(0) / kTau;

  res_.resize(0);
  for (unsigned int i = i_fft_; i < f_fft_ + 1; ++i){
    res_.push_back(psd_[i] - func.Eval(fftfreq_[i]));
  }

  return;
}

double Fid::CalcZeroCountFreq()
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
  //  double thresh = hyst_thresh_ * max;
  //  thresh = 10 * noise_;

  int i_zero = -1;
  int f_zero = -1;

  // iterate over vector
  for (unsigned int i = i_wf_; i < f_wf_; i++){

    // hysteresis check
    if (hyst){
      hyst = std::abs(wf_[i]) > hyst_thresh_ * env_[i];
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

  freq_ = 0.5 * (nzeros - 1) / (tf - ti);
  // todo: Fix this into a better error estimate. 
  freq_err_ = freq_ * sqrt(2) * (tm_[1] - tm_[0]) / (tf - ti);

  return freq_;
}

double Fid::CalcCentroidFreq()
{
  // Find the peak power
  double thresh = *std::max_element(psd_.begin()+1, psd_.end());
  thresh *= centroid_thresh_;

  // Find the indices for a window around the max
  int it_i = std::distance(psd_.begin(), 
    std::find_if(psd_.begin()+1, psd_.end(), 
      [thresh](double x) {return x > thresh;}));

  // reverse the iterators
  int it_f = -1 * std::distance(psd_.rend(),
    std::find_if(psd_.rbegin(), psd_.rend()-1, 
      [thresh](double x) {return x > thresh;}));

  // Now compute the power weighted average
  double pwfreq = 0.0;
  double pwfreq2 = 0.0;
  double pwsum = 0.0;

  for (int i = it_i; i < it_f; i++){
    pwfreq += psd_[i] * fftfreq_[i];
    pwfreq2 += psd_[i] * psd_[i] * fftfreq_[i];
    pwsum  += psd_[i];
  }

  freq_err_ = sqrt(pwfreq2 / pwsum - pow(pwfreq / pwsum, 2.0));
  freq_ = pwfreq / pwsum;

  return freq_;
}


double Fid::CalcAnalyticalFreq()
{
  // Set the fit function:
  // p0 -> fid frequency
  // p1 -> decay gamma
  // p2 -> amplitude
  // p3 -> baseline offset
  // p4 -> phi
  std::string fcn("[3]  + [2] * (([0] * cos([4]) + [1] * sin([4]))^2 ");
  fcn += std::string("+ (x * sin([1]))^2) / (([1]^2 + [0]^2 - x^2)^2 ");
  fcn += std::string("+ 4 * x^2 * [1]^2)");

  f_fit_ = TF1("f_fit_", fcn.c_str(), fftfreq_[i_fft_], fftfreq_[f_fft_]);

  // Set the parameter guesses
  for (unsigned int i = 0; i < 5; i++){
    f_fit_.SetParameter(i, guess_[i]);
  }

  // Limits
  f_fit_.SetParLimits(4, 0.0, kTau);

  FreqFit(f_fit_);

  freq_ = f_fit_.GetParameter(0);

  return freq_;
}


double Fid::CalcLorentzianFreq()
{
  // Set the fit function
  std::string fcn("[2] / (1 + ((x - [0]) / (0.5 * [1]))^2) + [3]");
  f_fit_ = TF1("f_fit_", fcn.c_str(), fftfreq_[i_fft_], fftfreq_[f_fft_]);

  // Set the parameter guesses
  for (int i = 0; i < 4; i++){
    f_fit_.SetParameter(i, guess_[i]);
  }

  FreqFit(f_fit_);

  freq_ = f_fit_.GetParameter(0);

  return freq_;
}


double Fid::CalcSoftLorentzianFreq()
{
  // Set the fit function
  f_fit_ = TF1("f_fit_", "[2] / (1 + ((x - [0]) / (0.5 * [1]))^[4]) + [3]");

  // Set the parameter guesses
  for (int i = 0; i < 5; i++){
    f_fit_.SetParameter(i, guess_[i]);
  }

  f_fit_.SetParLimits(4, 1.0, 3.0);

  FreqFit(f_fit_);

  freq_ = f_fit_.GetParameter(0);

  return freq_;
}


double Fid::CalcExponentialFreq()
{
  // Set the fit function
  f_fit_ = TF1("f_fit_", "[2] * exp(-abs(x - [0]) / [1]) + [3]");

  // Set the parameter guesses
  for (int i = 0; i < 4; i++){
    f_fit_.SetParameter(i, guess_[i]);
  }

  f_fit_.SetParameter(1, 0.5 * guess_[1] / std::log(2));

  FreqFit(f_fit_);

  freq_ = f_fit_.GetParameter(0);

  return freq_;
}


double Fid::CalcPhaseFreq(int poln)
{
  gr_time_series_ = TGraph(f_wf_ - i_wf_, &tm_[i_wf_], &phi_[i_wf_]);

  // Now set up the polynomial phase fit
  char fcn[20];
  sprintf(fcn, "pol%d", poln);
  f_fit_ = TF1("f_fit_", fcn);

  // Set the parameter guesses
  f_fit_.SetParameter(1, guess_[0] * kTau);

  // Adjust to ignore the edges
  int i = i_wf_ + edge_ignore_;
  int f = f_wf_ - edge_ignore_;

  // Do the fit.
  gr_time_series_.Fit(&f_fit_, "QMRSEX0", "C", tm_[i], tm_[f]);

  res_.resize(0);
  for (unsigned int idx = i; idx <= f; ++idx){
    res_.push_back(phi_[idx] - f_fit_.Eval(tm_[idx]));
  }

  freq_ = f_fit_.GetParameter(1) / kTau;
  freq_err_ = f_fit_.GetParError(1) / kTau;
  chi2_ = f_fit_.GetChisquare();

  return freq_;
}

double Fid::CalcPhaseDerivFreq(int poln)
{
  // Do the normal polynomial phase fit to set the fit function
  CalcPhaseFreq(poln);

  // Find the initial phase by looking at the function's derivative
  freq_ = f_fit_.Derivative(tm_[i_wf_]) / kTau;

  return freq_;
} 


double Fid::CalcSinusoidFreq()
{
  // Normalize the waveform by the envelope
  temp_.resize(wf_.size());

  std::transform(wf_.begin(), wf_.end(), env_.begin(), temp_.begin(),
    [](double x_wf, double x_env) {return x_wf / x_env;});

  gr_time_series_ = TGraph(f_wf_ - i_wf_, &tm_[i_wf_], &temp_[i_wf_]);    

  f_fit_ = TF1("f_fit_", "[1] * sin([0] * x) + [2] * cos([0] * x)");

  // Guess parameters
  f_fit_.SetParameter(0, kTau * CalcZeroCountFreq());
  f_fit_.SetParameter(1, 0.5);
  f_fit_.SetParameter(2, 0.5);
  f_fit_.SetParLimits(1, 0.0, 1.1);
  f_fit_.SetParLimits(2, 0.0, 1.1);

  // Adjust to ignore the edges
  int i = i_wf_ + edge_ignore_;
  int f = f_wf_ - edge_ignore_;

  // Do the fit.
  gr_time_series_.Fit(&f_fit_, "QMRSEX0", "C", tm_[i], tm_[f]);

  res_.resize(0);
  for (unsigned int i = i_fft_; i < f_fft_ + 1; ++i){
    res_.push_back(wf_[i] - f_fit_.Eval(tm_[i]));
  }

  freq_ = f_fit_.GetParameter(0) / kTau;
  freq_err_ = f_fit_.GetParError(0) / kTau;
  chi2_ = f_fit_.GetChisquare();

  return freq_;
}


Method Fid::ParseMethod(const std::string& m)
{
  using std::string;

  // Test each case iteratively, Zero count first.
  string str1("ZEROCOUNT");
  string str2("ZC");

  if (boost::iequals(m, str1) || boost::iequals(m , str2)) {
    return Method::ZC;
  }

  str1 = string("CENTROID");
  str2 = string("CN");

  if (boost::iequals(m, str1) || boost::iequals(m , str2)) {
    return Method::CN;
  }

  str1 = string("ANALYTICAL");
  str2 = string("AN");

  if (boost::iequals(m, str1) || boost::iequals(m , str2)) {
    return Method::AN;
  }

  str1 = string("LORENTZIAN");
  str2 = string("LZ");

  if (boost::iequals(m, str1) || boost::iequals(m , str2)) {
    return Method::LZ;
  }

  str1 = string("EXPONENTIAL");
  str2 = string("EX");

  if (boost::iequals(m, str1) || boost::iequals(m , str2)) {
    return Method::EX;
  }

  str1 = string("PHASE");
  str2 = string("PH");

  if (boost::iequals(m, str1) || boost::iequals(m , str2)) {
    return Method::PH;
  }

  str1 = string("SINUSOID");
  str2 = string("SN");

  if (boost::iequals(m, str1) || boost::iequals(m , str2)) {
    return Method::SN;
  }

  // If the method hasn't matched yet, use the current method and warn them.
  std::cout << "Warning: Method string not matched. " << std::endl;
  std::cout << "Method not changed." << std::endl;

  return freq_method_;
}


// Analyze an FID with all methods and print to output stream csv style.
void Fid::WriteFreqCsv(std::ofstream& out)
{
  // Test all the frequency extraction methods and write the results
  out << CalcZeroCountFreq() << ", " << 0.0 << ", ";
  out << CalcCentroidFreq() << ", " << 0.0 << ", ";
  out << CalcAnalyticalFreq() << ", " << chi2_ << ", ";
  out << CalcLorentzianFreq() << ", " << chi2_ << ", ";
  out << CalcSoftLorentzianFreq() << ", " << chi2_ << ", ";
  out << CalcExponentialFreq() << ", " << chi2_ << ", ";
  out << CalcPhaseFreq() << ", " << chi2_ << ", ";
  out << CalcSinusoidFreq() << ", " << chi2_ << std::endl;
}

void Fid::WriteMethodCsv(std::ofstream& out)
{
  using namespace std::chrono;

  // Test all the frequency extraction methods and
  // write the freq, ferr, chi2 and calculation time.

  long dt;
  auto t0 = high_resolution_clock::now();

  t0 = high_resolution_clock::now();
  GetFreq("ZC");
  dt = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  out << freq_ << ", " << freq_err_ << ", " << chi2_ << ", " << dt << ", ";

  t0 = high_resolution_clock::now();
  GetFreq("CN");
  dt = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  out << freq_ << ", " << freq_err_ << ", " << chi2_ << ", " << dt << ", ";

  t0 = high_resolution_clock::now();
  GetFreq("AN");
  dt = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  out << freq_ << ", " << freq_err_ << ", " << chi2_ << ", " << dt << ", ";

  t0 = high_resolution_clock::now();
  GetFreq("LZ");
  dt = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  out << freq_ << ", " << freq_err_ << ", " << chi2_ << ", " << dt << ", ";

  t0 = high_resolution_clock::now();
  GetFreq("EX");
  dt = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  out << freq_ << ", " << freq_err_ << ", " << chi2_ << ", " << dt << ", ";

  t0 = high_resolution_clock::now();
  GetFreq("PH");
  dt = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  out << freq_ << ", " << freq_err_ << ", " << chi2_ << ", " << dt << ", ";

  t0 = high_resolution_clock::now();
  GetFreq("SN");
  dt = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  out << freq_ << ", " << freq_err_ << ", " << chi2_ << ", " << dt;
  out << std::endl;

}

// Test all the frequency extraction methods using the phase
void Fid::WritePhaseFreqCsv(std::ofstream& out)
{
  out << CalcPhaseFreq() << ", " << chi2_ << ", ";
  out << CalcPhaseFreq(2) << ", " << chi2_ << ", ";
  out << CalcPhaseFreq(3) << ", " << chi2_ << ", ";
  out << CalcPhaseDerivFreq() << ", " << chi2_ << ", ";
  out << CalcPhaseDerivFreq(2) << ", " << chi2_ << ", ";
  out << CalcPhaseDerivFreq(3) << ", " << chi2_ << ", ";
  out << CalcSinusoidFreq() << ", " << chi2_ << std::endl;
}

} // fid
