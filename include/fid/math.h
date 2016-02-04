#ifndef FID_INCLUDE_FID_MATH_H_
#define FID_INCLUDE_FID_MATH_H_

/*===========================================================================*\

author: Matthias W. Smith
email: mwmsith2@uw.edu

notes: 

	This header defines some commonly used mathematical and digital 
	signal processing functions for the FID library.

\*===========================================================================*/

//--- std includes ----------------------------------------------------------//
#include <vector>
#include <functional>
#include <numeric>
#include <random>
#include <complex>
#include <ctgmath>
#include <cmath>

//--- other includes --------------------------------------------------------//
#include <armadillo>
#include "fftw3.h"

//--- project includes ------------------------------------------------------//
#include "fid/params.h"

namespace fid
{
// An alias
typedef std::complex<double> cdouble;

//--- linalg template functions ---------------------------------------------//

// A template function to handle vector addition.
template <typename T>
inline std::vector<T>& operator+(std::vector<T>& a, std::vector<T>& b)
{
  assert(a.size() == b.size());

  transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<T>());
  return a;
}

// A template function to handle vector subtraction.
template <typename T>
inline std::vector<T>& operator-(std::vector<T>& a, std::vector<T>& b)
{
  assert(a.size() == b.size());

  transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<T>());
  return a;
}

// A template function to handle vector multiplying with a scalar.
template <typename T>
inline std::vector<T>& operator*(T c, std::vector<T>& a)
{
  for (auto it = a.begin(); it != a.end(); ++it){
    *it = c * (*it);
  }

  return a;
}

template <typename T>
inline void floor(std::vector<T>& v) 
{
  for (auto it = v.begin(); it != v.end(); ++it) {
    *it = std::floor(*it);
  }
}

template <typename T>
inline void cross(const std::vector<T>& u, 
                  const std::vector<T>& v, 
                  std::vector<T>& res)
{
    res[0] = u[1] * v[2] - u[2] * v[1];
    res[1] = u[2] * v[0] - u[0] * v[2];
    res[2] = u[0] * v[1] - u[1] * v[0];
}

// Standard deviation calculation based on whole vector.
template <typename T>
inline double stdev(const T& v) {
    auto x1 = std::accumulate(v.begin(), v.end(), 0.0, 
        [](double x, double y) { return x + y; });

	auto x2 = std::accumulate(v.begin(), v.end(), 0.0, 
		[](double x, double y) { return x + y*y; });

  double N = (double) std::distance(v.begin(), v.end());
	return std::sqrt(x2/N - (x1/N) * (x1/N));
}

// Standard deviation calculation based on start/stop iterators.
template <typename T>
inline double stdev(const T& begin, const T& end) {
    auto x1 = std::accumulate(begin, end, 0.0, 
        [](double x, double y) { return x + y; });

	auto x2 = std::accumulate(begin, end, 0.0, 
		[](double x, double y) { return x + y*y; });

  double N = (double) std::distance(begin, end);
  return std::sqrt(x2/N - (x1/N) * (x1/N));
}

// Add white noise to an array.
template <typename T>
void addwhitenoise(std::vector<T>& v, T snr) {
  static std::default_random_engine gen(clock());
  std::normal_distribution<T> nrm;

  T mean = 0.0;
  T min = 0.0;
  T max = 0.0;
  for (auto it = v.begin(); it != v.end(); ++it) {

    if (*it > max) { 
      max = *it;
    } else {
      min = *it;
    }
    mean += *it;
  }

  // normalize to mean
  mean /= v.size();
  max -= mean;
  min = std::abs(min - mean);

  T amp = max > min ? max : min;
  T scale = amp / sqrt(snr);

  for (auto &x : v){
    x += nrm(gen) * scale;
  }
}

// Construct a range from params first, step, last
template<typename T>
inline std::vector<T> construct_range(const T& i, const T& f, const T& d) {
	std::vector<T> res;
	for (T x = i; x <= f; x += d){
	  res.push_back(x);
	}
	return res;
}

// Construct a range from vector <first, last, step>
template<typename T>
inline std::vector<T> construct_range(const std::vector<T> &range_vec) {
	std::vector<T> res;
	for (T x = range_vec[0]; x <= range_vec[1]; x += range_vec[2]){
	  res.push_back(x);
	}
	return res;
}

template<typename T>
inline std::vector<T> construct_linspace(const T& i, const T& f, const int& n=0) {
    std::vector<T> res;
    int N;
    if (n == 0) {
      N = f - i;
    } else {
      N = n;
    }
    double d = (f - i) / (n - 1);
    for (int j = 0; j < n; ++j) {
        res.push_back(i + j*d);
    }
    return res;
}

template<typename T>
inline std::vector<T> construct_linspace(const std::vector<T>& vals) {
    std::vector<T> res;
    int n = (int)(vals[2] + 0.5);
    double d = (vals[1] - vals[0]) / (n - 1);
    for (int j = 0; j < n; ++j) {
        res.push_back(vals[0] + j*d);
    }
    return res;
}

// Get a linear gradient with a max(abs(grad)) == 1 and mean(grad) == 0
std::vector<double> normalized_gradient(int npoints, int poln=1); 

namespace dsp
{
std::vector<cdouble> 
fft(const std::vector<cdouble>& v);

std::vector<cdouble> 
ifft(const std::vector<cdouble>& v);

std::vector<cdouble> rfft(const std::vector<double>& v);
std::vector<double> irfft(const std::vector<cdouble>& v);

std::vector<double> hilbert(const std::vector<double>& v);
std::vector<double> psd(const std::vector<double>& v);

std::vector<double> norm(const std::vector<double>& v);
std::vector<double> norm(const std::vector<cdouble>& v);

std::vector<double> fftfreq(const std::vector<double>& tm);
std::vector<double> fftfreq(const int N, const double dt);

std::vector<double> phase(const std::vector<double>& v);
std::vector<double> phase(const std::vector<double>& wf_re, 
                          const std::vector<double>& wf_im);

std::vector<double> envelope(const std::vector<double>& v);
std::vector<double> envelope(const std::vector<double>& wf_re, 
                             const std::vector<double>& wf_im);	


arma::mat wvd(const std::vector<double>& wf, bool upsample=false);
arma::cx_mat wvd_cx(const std::vector<double>& wf, bool upsample=false, const int window=0);
arma::cx_vec acorrelation(const arma::cx_vec &v, const int idx,  const int window=0);

std::vector<double> savgol3(const std::vector<double>& v);
std::vector<double> savgol5(const std::vector<double>& v);
std::vector<double> convolve(const std::vector<double>& v, 
                             const std::vector<double>& filter);
int convolve(const std::vector<double>& v, 
             const std::vector<double>& filter, 
             std::vector<double>& res);

template <typename T>
std::vector<T> lowpass(const std::vector<T>& v, double cut_idx=-1, int n=3) {
	// A simple Butterworth n-order filter.
	if (cut_idx == -1) cut_idx = v.size() / 2;
	std::vector<T> filtered_wf = v;

	std::transform(v.begin(), v.end(), filtered_wf.begin(), // lambda filter
				  [cut_idx, n](double x) { 
				  	return sqrt(1.0 / (1.0 + pow(x / cut_idx, 2*n))); 
				  });

	return filtered_wf;
}

template <typename T>
arma::Col<T> rconvolve(const arma::Col<T>& v, int idx=0) {
	int ridx = v.n_elem - idx;
	arma::Col<T> rv(arma::flipud(v));
	arma::Col<T> res(v.n_elem, arma::fill::zeros);

	if (idx > ridx) {
		std::transform(v.begin() + idx, v.end(), rv.begin(), res.begin(),
			[](T z1, T z2) { return z1 * std::conj(z2); });
	} else {
		std::transform(rv.begin() + ridx, rv.end(), v.begin(), res.begin(),
			[](T z2, T z1) { return z1 * std::conj(z2); });
	}
	return res;
} 

} // ::dsp

} // ::fid

#endif
