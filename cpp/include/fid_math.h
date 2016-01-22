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
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath> 
using std::vector;
using std::cout;
using std::endl;

//--- other includes --------------------------------------------------------//
#include <armadillo>
#include "fftw3.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"

namespace fid
{
//--- linalg template functions ---------------------------------------------//

// A template function to handle vector addition.
template <typename T>
inline vector<T>& operator+(vector<T>& a, vector<T>& b)
{
  assert(a.size() == b.size());

  transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<T>());
  return a;
}

// A template function to handle vector subtraction.
template <typename T>
inline vector<T>& operator-(vector<T>& a, vector<T>& b)
{
  assert(a.size() == b.size());

  transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<T>());
  return a;
}

// A template function to handle vector multiplying with a scalar.
template <typename T>
inline vector<T>& operator*(T c, vector<T>& a)
{
  for (auto it = a.begin(); it != a.end(); ++it){
    *it = c * (*it);
  }

  return a;
}

inline void cross(const vec& u, const vec& v, vec& res)
{
    res[0] = u[1] * v[2] - u[2] * v[1];
    res[1] = u[2] * v[0] - u[0] * v[2];
    res[2] = u[0] * v[1] - u[1] * v[0];
}

// Standard deviation calculation based on whole vector.
template <typename T>
inline double stdev(const T& wf) {
    auto x1 = std::accumulate(wf.begin(), wf.end(), 0.0, 
        [](double x, double y) { return x + y; });

	auto x2 = std::accumulate(wf.begin(), wf.end(), 0.0, 
		[](double x, double y) { return x + y*y; });

  double N = (double) std::distance(wf.begin(), wf.end());
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
void addnoise(vector<T>& wf, T snr, int seed=0) {
  static std::default_random_engine gen(seed);
  std::normal_distribution<T> nrm;

  T max = *std::max_element(wf.begin(), wf.end());
  T min = *std::min_element(wf.begin(), wf.end());
  T scale = max > min ? max : min;

  for (auto &x : wf){
    x += scale * nrm(gen) / sqrt(snr);
  }
}

// Construct a range from params first, step, last
template<typename T>
inline vector<T> construct_range(const T& i, const T& f, const T& d) {
	vector<T> res;
	for (T x = i; x <= f; x += d){
	  res.push_back(x);
	}
	return res;
}

// Construct a range from vector <first, last, step>
template<typename T>
inline vector<T> construct_range(const vector<T> &range_vec) {
	vector<T> res;
	for (T x = range_vec[0]; x <= range_vec[1]; x += range_vec[2]){
	  res.push_back(x);
	}
	return res;
}

template<typename T>
inline vector<T> construct_linspace(const T& i, const T& f, const int& n) {
    vector<T> res;
    double d = (f - i) / (n - 1);
    for (int j = 0; j < n; ++j) {
        res.push_back(i + j*d);
    }
    return res;
}

template<typename T>
inline vector<T> construct_linspace(const vector<T>& vals) {
    vector<T> res;
    int n = (int)(vals[2] + 0.5);
    double d = (vals[1] - vals[0]) / (n - 1);
    for (int j = 0; j < n; ++j) {
        res.push_back(vals[0] + j*d);
    }
    return res;
}

namespace dsp
{
cvec fft(const vec& wf);
vec ifft(const cvec& fft_vec);

vec hilbert(const vec& wf);
vec hilbert(cvec fft_vec);

vec psd(const vec& wf);
vec psd(const cvec& fft);

vec fftfreq(const vec& tm);
vec fftfreq(const int N, const double dt);

vec phase(const vec& wf);
vec phase(const vec& wf_re, const vec& wf_im);

vec envelope(const vec& wf);
vec envelope(const vec& wf_re, const vec& wf_im);	

arma::cx_mat wvd_cx(const vec& wf, bool upsample=false, const int window=0);
arma::mat wvd(const vec& wf, bool upsample=false);
arma::cx_vec acorrelation(const arma::cx_vec &v, const int idx,  const int window);

vec savgol3(const vec& wf);
vec savgol5(const vec& wf);
vec convolve(const vec& wf, const vec& filter);
int convolve(const vec& wf, const vec& filter, vec& res);

template <typename T>
vector<T> lowpass(const vector<T>& wf, double cut_idx=-1, int n=3) {
	// A simple Butterworth n-order filter.
	if (cut_idx == -1) cut_idx = wf.size() / 2;
	vector<T> filtered_wf = wf;

	std::transform(wf.begin(), wf.end(), filtered_wf.begin(), // lambda filter
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
