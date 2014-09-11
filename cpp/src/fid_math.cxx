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
  auto wf_vec = wf; // copy waveform since fftw destroys it.

  // Plan and execute the fft.
  fftw_plan wf_to_fft;
  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);
  wf_to_fft = fftw_plan_dft_r2c_1d(N, &wf_vec[0], fft_ptr, FFTW_ESTIMATE);
  fftw_execute(wf_to_fft);
  fftw_destroy_plan(wf_to_fft);

  for (auto& val: fft_vec) {
  	val /= Nroot;
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
  vec ifft_vec(n, 0.0);
  auto fft_vec = fft;

  // Plan and execute the fft.
  fftw_plan fft_to_wf;
  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);
  fft_to_wf = fftw_plan_dft_c2r_1d(N, fft_ptr, &ifft_vec[0], FFTW_ESTIMATE);
  fftw_execute(fft_to_wf);
  fftw_destroy_plan(fft_to_wf);

  // fftw is unnormalized, so we need to fix that.
  for (auto& val: ifft_vec) {
  	val /= Nroot;
  }

  return ifft_vec;
}

vec dsp::hilbert(const vec& wf)
{
	// Return the call to the fft version.
	auto fft_vec = dsp::fft(wf);
	return dsp::hilbert(fft_vec);
}

vec dsp::hilbert(cvec& fft_vec)
{
	// Multiply in the -i.
	for (auto& val : fft_vec) {
		val = std::complex<double>(-val.imag(), val.real());
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

} // ::fid