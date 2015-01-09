#include "fid_math.h"

namespace fid 
{

cvec dsp::fft(const vec &wf)
{
  // Grab some useful constants.
  int N = wf.size();  
  int n = N / 2 + 1;  // size of rfft
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  cvec fft_vec(n, 0.0);
  auto wf_vec = wf; // copy waveform since fftw destroys it

  // Plan and execute the fft.
  fftw_plan wf_to_fft;
  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);
  wf_to_fft = fftw_plan_dft_r2c_1d(N, &wf_vec[0], fft_ptr, FFTW_ESTIMATE);
  fftw_execute(wf_to_fft);
  fftw_destroy_plan(wf_to_fft);

  for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
    *it /= Nroot;
  }

  return fft_vec;
}

vec dsp::ifft(const cvec& fft)
{
  // Grab some useful constants.
  int n = fft.size();
  int N = 2 * (n - 1);
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  vec ifft_vec(N, 0.0);
  fftw_complex *fft_ptr = new fftw_complex[n];
  memcpy(fft_ptr, &fft[0], sizeof(fftw_complex) * n);

  // Plan and execute the fft.
  fftw_plan fft_to_wf;
  fft_to_wf = fftw_plan_dft_c2r_1d(N, fft_ptr, &ifft_vec[0], FFTW_ESTIMATE);
  fftw_execute(fft_to_wf);
  fftw_destroy_plan(fft_to_wf);

  // fftw is unnormalized, so we need to fix that.
  for (auto it = ifft_vec.begin(); it != ifft_vec.end(); ++it) {
  	*it /= Nroot;
  }

  return ifft_vec;
}

vec dsp::hilbert(const vec& wf)
{
	// Return the call to the fft version.
	auto fft_vec = dsp::fft(wf);
	return dsp::hilbert(fft_vec);
}

vec dsp::hilbert(cvec fft_vec)
{
	// Multiply in the -i.
	for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
		*it = std::complex<double>(-(*it).imag(), (*it).real());
	}

	// Reverse the fft.
	return ifft(fft_vec);
}

vec dsp::psd(const vec& wf)
{
	return dsp::psd(dsp::fft(wf));
}

vec dsp::psd(const cvec& fft_vec)
{
	// Instatiate the power vector and fill it with the magnitude of fft_vec.
	vec power(fft_vec.size(), 0.0);

	for (int i = 0; i < fft_vec.size(); ++i) {
		power[i] = std::norm(fft_vec[i]);
	}

	return power;
}

// Helper function to get frequencies for FFT
vec dsp::fftfreq(const vec& tm) 
{
	int N = tm.size();
	double dt = (tm[N-1] - tm[0]) / (N - 1); // sampling rate

	return dsp::fftfreq(N, dt);
}

vec dsp::fftfreq(const int N, const double dt)
{
	// Instantiate return vector.
	vec freq;

	// Handle both even and odd cases properly.
	if (N % 2 == 0) {

		freq.resize(N/2 + 1);
		
		for (int i = 0; i < N/2 + 1; ++i) {
			freq[i] = i / (dt * N);
		}

	} else {

		freq.resize((N + 1) / 2);

		for (int i = 0; i < (N + 1) / 2 + 1; ++i){
			freq[i] = i / (dt * N);
		}
	}

	return freq;
}

// Calculates the phase by assuming the real signal is harmonic.
vec dsp::phase(const vec& wf)
{
	return dsp::phase(wf, dsp::hilbert(wf));
}

vec dsp::phase(const vec& wf_re, const vec& wf_im)
{
	vec phase(wf_re.size(), 0.0);

	// Calculate the modulo-ed phase
  	std::transform(wf_re.begin(), wf_re.end(), wf_im.begin(), phase.begin(),
    			   [](double re, double im) { return std::atan2(im, re); });

	// Now unwrap the phase
	double thresh = params::max_phase_jump;
	int k = 0; // to track the winding number
  	for (auto it = phase.begin(); it != phase.end(); ++it) {

    	// Add current total
    	*it += k * kTau;

    	// Check for large jumps, both positive and negative.
    	if (*(it - 1) - *it > thresh) {
	      	k++;
    	  	*it += kTau;

	    } else if (*it - *(it - 1) > thresh) {
	        k--;
	        *it -= kTau;
	    }
    }

    return phase;
}

vec dsp::envelope(const vec& wf)
{
	return dsp::envelope(wf, dsp::hilbert(wf));
}

vec dsp::envelope(const vec& wf_re, const vec& wf_im)
{
 	// Set the envelope function
	vec env(wf_re.size(), 0.0);

	std::transform(wf_re.begin(), wf_re.end(), wf_im.begin(), env.begin(),
    			   [](double r, double i) { return std::sqrt(r*r + i*i); });

	return env;
}


arma::mat dsp::wvd(const vec& wf) 
{
  // Instiate the return matrix (and autoconvultion matrix-Ronaldo)
  arma::mat res(wf.size(), wf.size(), arma::fill::zeros);

  // Make the signal harmonic
  arma::cx_vec v(wf.size());

  auto wf_im = dsp::hilbert(wf);

  for (int i = 0; i < wf.size(); ++i) {
    v[i] = arma::cx_double(wf[i], wf_im[i]);
  }

  // Now compute the Wigner-Ville Distribution
  int idx = 0;
  for (auto it = wf.begin(); it != wf.end(); ++it) {
    arma::vec fft(arma::abs(arma::fft(dsp::rconvolve(v, idx))));
    res.col(idx++) = fft;
 
  }

  return res;
}

arma::cx_mat dsp::auto_correlation(const vec& wf)
{
  arma::mat real_ACF(wf.size(),wf.size(), arma::fill::zeros);
  arma::mat im_ACF(wf.size(), wf.size(), arma::fill::zeros);
  
  arma::cx_vec analytic_wf(wf.size());
  
  auto wf_im = dsp::hilbert(wf);

  for (int i=0; i<wf.size(); ++i) {
    analytic_wf[i] = arma::cx_double(wf[i], wf_im[i]);
  }

  int tao = 0;
  for (auto it = wf.begin(); it !=wf.end(); ++it) {
    arma::cx_vec auto_corr = dsp::rconvolve(analytic_wf, tao);
    // im_wvd.col(tao++) = imag(fft_ACF);
    // real_wvd.col(tao++) = real(fft_ACF);
    real_ACF.col(tao) = arma::real(auto_corr);
    im_ACF.col(tao++) = arma::imag(auto_corr);
  }
  arma::cx_mat ACF(real_ACF,im_ACF);
    
  return ACF;
    
}

arma::cx_mat dsp::complex_wvd(const vec& wf)
{
  arma::cx_mat cx_wvd(wf.size(), wf.size(), arma::fill::zeros);

  arma::cx_vec analytic_wf(wf.size());
  
  auto wf_im = dsp::hilbert(wf);

  for (int i=0; i<wf.size(); ++i) {
    analytic_wf[i] = arma::cx_double(wf[i], wf_im[i]);
  }

  int tao = 0;
  for (auto it = wf.begin(); it !=wf.end(); ++it) {
    // arma::cx_vec fft(arma::fft(dsp::rconvolve(analytic_wf, tao)));
    cx_wvd.col(tao++) = analytic_wf;
  }

  return cx_wvd;
}

std::pair<arma::cx_mat, arma::cx_mat> dsp::wvd_associates(const vec& wf)
{
  arma::cx_mat cx_wvd(2*wf.size(), wf.size(), arma::fill::zeros);
  arma::cx_mat ACF(2*wf.size(), wf.size(), arma::fill::zeros);
  
  vec wf_upsample(2*wf.size());
  arma::cx_vec v(4*wf.size());

  auto it1 = wf_upsample.begin();
  for (auto it2 = wf.begin(); it2 != wf.end(); ++it2) {
    *(it1++) = *it2;
    *(it1++) = *it2;
  }

  auto wf_im = dsp::hilbert(wf_upsample);
  
  for (uint i = 0; i < 4*wf.size(); ++i) {
    v[i] = arma::cx_double(wf_upsample[i], wf_im[i]);
    if (i >= 2*wf.size()){
      v[i] = arma::cx_double(wf_upsample[i-2*wf.size()], wf_im[i-2*wf.size()]);
    }
  }

  int tao = 0;
  arma::cx_vec fft_i(wf.size());
  for (auto it = wf.begin(); it !=wf.end(); ++it) {
    // arma::cx_vec auto_corr = dsp::rconvolve_periodic(v, tao);
    ACF.col(tao) = dsp::rconvolve_periodic(v,tao);
    cx_wvd.col(tao++)= arma::fft(ACF.col(tao), v.n_elem/2);
    /* fft_i.zeros();
    for (int k=0; k<wf.size(); k++){
      fft_i[k] =arma::cx_double(fft[k].real(), fft[k].imag());
    }
    cx_wvd.col(tao) = fft_i;*/
  }

  return std::make_pair(ACF, cx_wvd);
}

std::pair<arma::cx_mat, arma::cx_mat> dsp::wvd_dsr(const vec& wf)
{
   // Instiate the return matrix
  arma::cx_mat res(2*wf.size(), wf.size(), arma::fill::zeros);
  arma::cx_mat acf(2*wf.size(), wf.size(), arma::fill::zeros);

  // Artificially double the sampling rate by doubling each sample.
  vec wf_upsample(wf.size() * 2, 0.0);

  auto it1 = wf_upsample.begin();
  for (auto it2 = wf.begin(); it2 != wf.end(); ++it2) {
    *(it1++) = *it2;
    *(it1++) = *it2;
  }

  // Make the signal harmonic
  arma::cx_vec v(wf_upsample.size());

  auto wf_im = dsp::hilbert(wf_upsample);

  for (uint i = 0; i < wf_upsample.size(); ++i) {
    v[i] = arma::cx_double(wf_upsample[i], wf_im[i]);
  }

  // Now compute the Wigner-Ville Distribution
  for (int idx = 0; idx < v.n_elem/2; ++idx) {
    acf.col(idx) = dsp::rconvolve(v, idx);
    res.col(idx) = arma::fft(acf.col(idx));
  }

  return std::make_pair(acf, res);
}

} // ::fid
