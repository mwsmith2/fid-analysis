#include "fid/math.h"
#include "TMultiGraph.h"
#include "fid.h"

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
  int n = N / 2 + 1;  // size of rfft
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<cdouble> fft_vec(N, cdouble(0.0, 0.0));
  auto wfm_vec = v; // copy waveform since fftw destroys it
  //Note forward FFT's do not modify the input vector as long as its not
  //the same as the ouput vector

  // Plan and execute the fft.
  fftw_complex *fft_ptr = reinterpret_cast<fftw_complex *>(&fft_vec[0]);
  fftw_complex *wfm_ptr = reinterpret_cast<fftw_complex *>(&wfm_vec[0]);  

  fftw_plan wf_to_fft = fftw_plan_dft_1d(N, 
                                         wfm_ptr,
                                         fft_ptr,
                                         FFTW_FORWARD,
                                         FFTW_ESTIMATE);  

  fftw_execute(wf_to_fft);
  fftw_destroy_plan(wf_to_fft);

  for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
    *it /= Nroot;
  }
  
  return fft_vec;
}


std::vector<cdouble> dsp::ifft(const std::vector<cdouble>& v)
{
  // Grab some useful constants.
  int n = v.size();
  int N = 2 * (n - 1);
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<cdouble> ifft_vec(N, cdouble(0.0, 0.0));

  fftw_complex *fft_ptr = new fftw_complex[n];
  fftw_complex *wfm_ptr;

  memcpy(fft_ptr, &v[0], sizeof(fftw_complex) * n);
  wfm_ptr = reinterpret_cast<fftw_complex *>(&ifft_vec[0]);  

  // Plan and execute the fft.
  fftw_plan fft_to_wf = fftw_plan_dft_1d(N, 
                                         fft_ptr,
                                         wfm_ptr,
                                         FFTW_BACKWARD,
                                         FFTW_ESTIMATE);  

  fftw_execute(fft_to_wf);
  fftw_destroy_plan(fft_to_wf);

  // fftw is unnormalized, so we need to fix that.
  for (auto it = ifft_vec.begin(); it != ifft_vec.end(); ++it) {
    *it /= Nroot;
  }

  delete[] fft_ptr;

  return ifft_vec;
}
std::vector<cdouble> dsp::rfft(const std::vector<double> &v)
{
  // Grab some useful constants.
  int N = v.size();  
  int n = N / 2 + 1;  // size of rfft
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<cdouble> fft_vec(n, 0.0);
  auto wf_vec = v; // copy waveform since fftw destroys it

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


std::vector<double> dsp::irfft(const std::vector<cdouble>& fft)
{
  // Grab some useful constants.
  int n = fft.size();
  int N = 2 * (n - 1);
  double Nroot = std::sqrt(N);

  // Instantiate the result vector.
  std::vector<double> ifft_vec(N, 0.0);
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

  delete[] fft_ptr;

  return ifft_vec;
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

 std::vector<double> dsp::hilbert_rot(const std::vector<double>&v ,const cdouble i)
{
  //Return the call to the fft version.
  auto fft_vec = dsp::rfft(v);

  //Multiply by i^2.
  for (auto it = fft_vec.begin(); it != fft_vec.end(); ++it) {
    if (it ==fft_vec.begin()) *it = cdouble(0,0);
    *it = i*(*it);//cdouble((*it).real(), (*it).imag());
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
  double thresh = params::max_phase_jump;
  double m = 0.0;

  int k = 0; // to track the winding number
  for (auto it = phase.begin() + 1; it != phase.end(); ++it) {
    
    // Add current total
    *it += k * kTau;
    m = *(it) - *(it - 1);
    
    // Check for large jumps, both positive and negative.
    while (std::abs(m) > thresh) {
      
      if (-m > thresh) {
        
        if (m + kTau > thresh) {
          std::cout << "Warning: jump over threshold." << std::endl;
          break;
        }
        
        k++;
        *it += kTau;
        
      } else if (m > thresh) {
        
        if (m - kTau < -thresh) {
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

  arma::cx_vec dsp::acorrelation(const arma::cx_vec &v, int idx , int window)
{
  int N = v.size();
  std::cout<< "Size of v is "<< N << std::endl;

  //Define Window function, with odd length to ensure even/odd symmetry
  arma::vec  win_func(window);
  arma::cx_vec conj_v = arma::conj(v);
  win_func.ones();

  arma::cx_vec acf(N);
  acf.zeros();
  std::cout << "size of acf vector is "<<acf.size()<<std::endl;

  int taumax = std::min({idx, std::abs(N/2-idx), window});
  //make sure taumax is odd
  if ((taumax % 2 ==0) &(taumax !=0)) taumax -= 1;
  std::cout<<"correlation range is " <<taumax<< std::endl;

  //@Todo: tune acf to give strictly real result and symmetric fft peak. 
  for (int i =-taumax+1; i < (taumax) ; i++) {
    //   std::cout<< "correlation variable is " <<i <<" idx is "<< idx<<std::endl;
    //  if (i == 0) acf[taumax] = .5*(conj_v(idx-i)*v(idx+i)+conj_v(idx+i)*v(idx-i));
    if ((i < 0)&(idx+taumax<N)) {
      acf(N+i) = conj_v(idx-i)*v(idx+i);//starts filling acf at index 0
    }   
    if ((i  >= 0 )&(idx+taumax<N)) {
      acf(i) = conj_v(idx-i)*v(idx+i);
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

  return acf;
} 

std::vector<cdouble> dsp::wvd_prep(const std::vector<double>& wf, bool upsample, const int window)
{
  //Artificially double sampling rate if desired
   int M, N;
  if (upsample) {

    M = 2 * wf.size();
    N = wf.size();

  } else {

    M = wf.size();
    N = wf.size();
  }
  std::vector<double> wf_re(M, 0.0);

  auto it1 = wf_re.begin();
  for (auto it2 = wf.begin(); it2 != wf.end(); ++it2) {
    *(it1++) = *it2;
    if (upsample) {
      *(it1++) = *it2;
    }
  }

  //Make signal analytic and centered about zero.
  std::vector<double> wf_cen(M, 0.0);
  cdouble r = cdouble(1.0, 0.0);
  std::vector<double> wf_pi = dsp::hilbert_rot(wf_re,r);
  std::transform( wf_re.begin(), wf_re.end(), wf_pi.begin(), wf_cen.begin(),
                  [](double re, double pi) {return re - std::abs(re-pi);});

  auto wf_im = dsp::hilbert(wf_cen);
  std::vector<cdouble> v(M, cdouble(0.0,0.0));

  for (int i = 0; i < wf_re.size(); ++i) {
    v[i] = cdouble(wf_cen[i], wf_im[i]);
    //  phase[i] = (1.0 * i) / M * M_PI;
  }

  return v;
}

std::vector<double> dsp::WvdFreqExt(const std::vector<double>& wf,  bool upsample, const int window)
{
  int M, N;
  if (upsample) {

    M = 2 * wf.size();
    N = wf.size();

  } else {

    M = wf.size();
    N = wf.size();
  }

  // Initiate the return matrix
  arma::vec re_fft_acorr(N);
  std::vector<double> wvd(N, 0.0);
  

  // Artificially double the sampling rate by repeating each sample.
  std::vector<double> wf_re(M, 0.0);

  auto it1 = wf_re.begin();
  for (auto it2 = wf.begin(); it2 != wf.end(); ++it2) {
    *(it1++) = *it2;
    if (upsample) {
      *(it1++) = *it2;
    }
  }

  // Make the signal harmonic
  arma::cx_vec v(M);
  arma::vec phase(M);

  auto wf_im = dsp::hilbert(wf_re);
  // auto wf_re_cen = dsp::hilbert(wf_im);

  for (uint i = 0; i < M; ++i) {
    v[i] = arma::cx_double(wf_re[i], wf_im[i]);
    phase[i] = (1.0 * i) / M * M_PI;
  }
  
  //Get Frequency dimensionality for FFT
  //@ToDo give wvd proper frequency dimensionality functionability.
  std::vector<double> freq= dsp::fftfreq(N, .001);

  // Now compute the Wigner-Ville Distribution
  for (int idx =  window; idx < N-window; ++idx) {
    re_fft_acorr =  arma::real(arma::fft(dsp::acorrelation(v, idx, window)));
    arma::uword max_bin;
    double max_val = re_fft_acorr.max(max_bin);

    int i_range=0; int f_range = 0;
    int n = re_fft_acorr.size();
    for (int i =1 ; i<n+1; i++) {
      if (re_fft_acorr[max_bin+i]>.15*re_fft_acorr[max_bin]) {
        f_range= max_bin + i+1;
        n = f_range-i_range;
      }
      if ((max_bin-i>0)&(re_fft_acorr[max_bin-i]>.15*re_fft_acorr[max_bin])) {
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
      sum += .5*(re_fft_acorr[i+1]+re_fft_acorr[i])*dx;
    }
  
    for (int i = 0; i<peak.size(); i++) {
      peak[i]=re_fft_acorr[i_range+i]/sum;
    }
    
    //Now perform fit to Gaussian
    std::string gaussian("[2]*exp(-(x-[0])^2/(2*[1]^2))");
    TGraph gr3 = TGraph(range, &freq[i_range], &peak[0]);
    TF1 fit_func = TF1("fit_func", gaussian.c_str(), freq[i_range],freq[f_range-1]);
    
    fit_func.SetParameter(0, freq[max_bin]);
    fit_func.SetParameter(1, 2);
    fit_func.SetParameter(2, peak.max());
    //   fit_func.SetParameter(3, 0);
    
    gr3.Fit(&fit_func, "R");

    std::cout<<"Center Frequency of FFT is: " << fit_func.GetParameter(0)/2<<std::endl;
    
    wvd[idx] = fit_func.GetParameter(0)/2;
  }

  return wvd;
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
