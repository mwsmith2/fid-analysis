#include "fid/math.h"

namespace fid {

std::vector<double> normalized_gradient(int npoints, int poln)
{
  std::vector<double> grad;

  // First get the spacing right.
  for (int i = 0; i < npoints; i++){
    grad.push_back(pow(i, poln));
  }

  // Subtract off the average.
  double avg = std::accumulate(grad.begin(), grad.end(), 0.0) / grad.size();

  for (uint i = 0; i < grad.size(); i++){
    grad[i] -= avg;
  }

  // Normalize by largest value.
  double max = *std::max_element(grad.begin(), grad.end());

  for (uint i = 0; i < grad.size(); i++){
    grad[i] /= max;
  }

  return grad;
}

std::vector<cdouble> dsp::fft(const std::vector<cdouble> &v)
{
  // Grab some useful constants.
  int N = v.size();  
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<cdouble> fft_vec(N, cdouble(0.0, 0.0));
  auto wfm_vec = v; // copy waveform since fftw destroys it

  // Plan and execute the fft (-1 == exponent).
  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);
  fftw_complex *wfm_ptr = reinterpret_cast<fftw_complex *>(&wfm_vec[0]);  

  auto plan = fftw_plan_dft_1d(N, wfm_ptr, fft_ptr, -1, FFTW_ESTIMATE);  
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
    *it /= Nroot;
  }

  return fft_vec;
}

std::vector<cdouble> dsp::ifft(const std::vector<cdouble>& v)
{
  // Grab some useful constants.
  int N = v.size();
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<cdouble> wfm_vec(N, cdouble(0.0, 0.0));
  auto fft_vec = v;

  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);  
  fftw_complex *wfm_ptr = reinterpret_cast<fftw_complex *>(&wfm_vec[0]);  

  // Plan and execute the inverse fft (+1 == exponent).
  auto plan = fftw_plan_dft_1d(N, fft_ptr, wfm_ptr, +1, FFTW_ESTIMATE);  
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // fftw is unnormalized, so we need to fix that.
  for (auto it = wfm_vec.begin(); it != wfm_vec.end(); ++it) {
    *it /= Nroot;
  }

  return wfm_vec;
}

std::vector<cdouble> dsp::rfft(const std::vector<double> &v)
{
  // Grab some useful constants.
  int N = v.size();  
  int n = N / 2 + 1;  // size of rfft
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<cdouble> fft_vec(n, 0.0);
  auto wfm_vec = v; // copy waveform since fftw destroys it

  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);

  // Plan and execute the fft.
  auto plan = fftw_plan_dft_r2c_1d(N, &wfm_vec[0], fft_ptr, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
    *it /= Nroot;
  }

  return fft_vec;
}


std::vector<double> dsp::irfft(const std::vector<cdouble>& fft)
{
  // Grab some useful constants.
  int n = fft.size();
  int N = 2 * (n - 1);
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<double> wfm_vec(N, 0.0);
  std::vector<cdouble> fft_vec = fft;

  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);

  // Plan and execute the fft.
  auto plan = fftw_plan_dft_c2r_1d(N, fft_ptr, &wfm_vec[0], FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  // fftw is unnormalized, so we need to fix that.
  for (auto it = wfm_vec.begin(); it != wfm_vec.end(); ++it) {
  	*it /= Nroot;
  }

  return wfm_vec;
}

std::vector<double> dsp::hilbert(const std::vector<double>& v)
{
	// Return the call to the fft version.
	auto fft_vec = dsp::rfft(v);

  // Multiply in the -i.
  for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
    *it = cdouble((*it).imag(), -(*it).real());
  }

  // Reverse the fft.
  return irfft(fft_vec);
}

std::vector<double> dsp::psd(const std::vector<double>& v)
{
  // Perform fft on the original data.
	auto fft_vec = dsp::rfft(v);

  // Get the norm of the fft as that is the power.
	return dsp::norm(fft_vec);
}

std::vector<double> dsp::norm(const std::vector<double>& v)
{
  // Allocate the memory
  std::vector<double> res;
  res.reserve(v.size());

  // Iterate and push back the norm.
  for (auto it = v.begin(); it < v.end(); ++it) {
    res.push_back(std::norm(*it));
  }

  return res;
}

std::vector<double> dsp::norm(const std::vector<cdouble>& v)
{
  // Allocate the memory
  std::vector<double> res;
  res.reserve(v.size());

  // Iterate and push back the norm.
  for (auto it = v.begin(); it < v.end(); ++it) {
    res.push_back(std::norm(*it));
  }

  return res;
}

// Helper function to get frequencies for FFT
std::vector<double> dsp::fftfreq(const std::vector<double>& tm) 
{
	int N = tm.size();
	double dt = (tm[N-1] - tm[0]) / (N - 1); // sampling rate

	return dsp::fftfreq(N, dt);
}

std::vector<double> dsp::fftfreq(const int N, const double dt)
{
	// Instantiate return vector.
	std::vector<double> freq;

	// Handle both even and odd cases properly.
	if (N % 2 == 0) {

		freq.resize(N/2 + 1);
		
		for (int i = 0; i < freq.size(); ++i) {
			freq[i] = i / (dt * N);
		}

	} else {

		freq.resize((N + 1) / 2);

		for (int i = 0; i < freq.size(); ++i){
			freq[i] = i / (dt * N);
		}
	}

	return freq;
}

// Calculates the phase by assuming the real signal is harmonic.
std::vector<double> dsp::phase(const std::vector<double>& v)
{
	return dsp::phase(v, dsp::hilbert(v));
}

std::vector<double> dsp::phase(const std::vector<double>& wf_re, 
                               const std::vector<double>& wf_im)
{
  std::vector<double> phase(wf_re.size(), 0.0);
  
  // Calculate the modulo-ed phase
  std::transform(wf_re.begin(), wf_re.end(), wf_im.begin(), phase.begin(),
                 [](double re, double im) { return std::atan2(im, re); });
  
  // Now unwrap the phase
  double m = 0.0;

  int k = 0; // to track the winding number
  for (auto it = phase.begin() + 1; it != phase.end(); ++it) {
    
    // Add current total
    *it += k * kTau;
    m = *(it) - *(it - 1);
    
    // Check for large jumps, both positive and negative.
    while (std::abs(m) > kMaxPhaseJump) {
      
      if (-m > kMaxPhaseJump) {
        
        if (m + kTau > kMaxPhaseJump) {
          std::cout << "Warning: jump over threshold." << std::endl;
          break;
        }
        
        k++;
        *it += kTau;
        
      } else if (m > kMaxPhaseJump) {
        
        if (m - kTau < -kMaxPhaseJump) {
          std::cout << "Warning: jump over threshold." << std::endl;
          break;
        }
        
        k--;
        *it -= kTau;
      }

      m = *(it) - *(it - 1);
    }
  }
  
  return phase;
}

std::vector<double> dsp::envelope(const std::vector<double>& v)
{
	return dsp::envelope(v, dsp::hilbert(v));
}

std::vector<double> dsp::envelope(const std::vector<double>& wf_re, const std::vector<double>& wf_im)
{
 	// Set the envelope function
	std::vector<double> env(wf_re.size(), 0.0);

	std::transform(wf_re.begin(), wf_re.end(), wf_im.begin(), env.begin(),
    			   [](double r, double i) { return std::sqrt(r*r + i*i); });

	return env;
}

arma::cx_mat dsp::wvd_cx(const std::vector<double>& v, bool upsample)
{
  int M, N;
  if (upsample) {

    M = 2 * v.size();
    N = v.size();

  } else {

    M = v.size();
    N = v.size();
  }

  // Initiate the return matrix
  arma::cx_mat res(M, N, arma::fill::zeros);

  // Artificially double the sampling rate by repeating each sample.
  std::vector<double> wf_re(M, 0.0);

  auto it1 = wf_re.begin();
  for (auto it2 = v.begin(); it2 != v.end(); ++it2) {
    *(it1++) = *it2;
    if (upsample) {
      *(it1++) = *it2;
    }
  }

  // Make the signal harmonic
  arma::cx_vec v2(M);
  arma::vec phase(M);

  auto wf_im = dsp::hilbert(wf_re);

  for (uint i = 0; i < M; ++i) {
    v2[i] = arma::cx_double(wf_re[i], wf_im[i]);
    phase[i] = (1.0 * i) / M * M_PI;
  }

  // Now compute the Wigner-Ville Distribution
  for (int idx = 0; idx < N; ++idx) {
    res.col(idx) = arma::fft(dsp::rconvolve(v2, idx));
  }

  return res;
}

arma::mat dsp::wvd(const std::vector<double>& v, bool upsample)
{
  int M, N;
  if (upsample) {

    M = 2 * v.size();
    N = v.size();

  } else {

    M = v.size();
    N = v.size();
  }

  // Instiate the return matrix
  arma::mat res(M, N, arma::fill::zeros);

  // Artificially double the sampling rate by repeating each sample.
  std::vector<double> wf_re(M, 0.0);

  auto it1 = wf_re.begin();
  for (auto it2 = v.begin(); it2 != v.end(); ++it2) {
    *(it1++) = *it2;
    if (upsample) {
      *(it1++) = *it2;
    }
  }

  // Make the signal harmonic
  arma::cx_vec v2(M);

  auto wf_im = dsp::hilbert(wf_re);

  for (int i = 0; i < M; ++i) {
    v2[i] = arma::cx_double(wf_re[i], wf_im[i]);
  }

  // Now compute the Wigner-Ville Distribution
  for (int idx = 0; idx < N; ++idx) {
    res.col(idx) = arma::real(arma::fft(dsp::rconvolve(v2, idx))) ;
  }

  return res;
}

std::vector<double> dsp::savgol3(const std::vector<double>& v)
{
  std::vector<double> res(0, v.size());
  std::vector<double> filter = {-2.0, 3.0, 6.0, 7.0, 6.0, 3.0, -2.0};
  filter = (1.0 / 21.0) * filter;

  if (dsp::convolve(v, filter, res) == 0) {

    return res;

  } else {

    return v;
  }
}

std::vector<double> dsp::savgol5(const std::vector<double>& v)
{
  std::vector<double> res(0, v.size());
  std::vector<double> filter = {15.0, -55.0, 30.0, 135.0, 179.0, 135.0, 30.0, -55.0, 15.0};
  filter = (1.0 / 429.0) * filter;

  if (dsp::convolve(v, filter, res) == 0) {

    return res;

  } else {

    return v;
  }
}

int dsp::convolve(const std::vector<double>& v, const std::vector<double>& filter, std::vector<double>& res)
{
  int k = filter.size();
  int N = v.size();
  res.resize(N);

  // Check to make sure we can do something.
  if (N < k) {
    return -1;
  }

  // First take care of the beginning and end.
  for (int i = 0; i < k + 1; ++i) {
    res[i] = 0.0;
    res[N -1 - i] = 0.0;

    for (int j = i; j < i + k; ++j) {

      res[i] += v[abs(j - k/2)] * filter[j - i];
      res[N - 1 - i] += v[N - 1 - abs(k/2 - j)] * filter[j - i];
    }
  }

  // Now the rest of the elements.
  for (auto it = v.begin(); it != v.end() - k; ++it) {
    double val = std::inner_product(it, it + k, filter.begin(), 0.0);
    res[std::distance(v.begin(), it + k/2)] = val;
  }

  return 0;
}

std::vector<double> dsp::convolve(const std::vector<double>& v, const std::vector<double>& filter)
{
  std::vector<double> res(0.0, v.size());

  if (dsp::convolve(v, filter, res) == 0) {

    return res;

  } else {

    return v;
  }
}
 
} // ::fid
