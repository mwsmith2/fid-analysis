#ifndef FID_INCLUDE_FID_MATH_H_
#define FID_INCLUDE_FID_MATH_H_

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
#include "fftw3.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"

// This header defines some commonly used mathematical and digital 
// signal processing functions for the FID library.

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

// Standard deviation calculation based on whole vector.
template <typename T>
T stdev(const T& wf) {
	auto x = std::accumulate(wf.begin(), wf.end(), 0.0, 
		[](double x, double y) { return x + y*y; });

	return std::sqrt(x);
}

// Standard deviation calculation based on start/stop iterators.
template <typename T>
T stdev(const T& begin, const T& end) {
	auto x = std::accumulate(begin, end, 0.0, 
		[](double x, double y) { return x + y*y; });

	return std::sqrt(x);
}

namespace dsp
{
cvec fft(const vec& wf);
vec ifft(const cvec& fft_vec);

vec hilbert(const vec& wf);
vec hilbert(cvec& fft_vec);

vec psd(const vec& wf);
vec psd(const cvec& fft);

vec fftfreq(const vec& tm);
vec fftfreq(const int N, const double dt);

vec phase(const vec& wf);
vec phase(const vec& wf_re, const vec& wf_im);

vec envelope(const vec& wf);
vec envelope(const vec& wf_re, const vec& wf_im);	

template <typename T>
std::vector<T> lowpass(const std::vector<T>& wf, double cut_idx=0, int n=3) {
	// A simple Butterworth n-order filter.
	std::vector<T> filtered_wf = wf;

	std::transform(wf.begin(), wf.end(), filtered_wf.begin(), // lambda filter
				  [cut_idx, n](double x) { 
				  	return sqrt(1.0 / (1.0 + pow(x / cut_idx, 2*n))); 
				  });

	return filtered_wf;
}

} // ::dsp

} // ::fid

#endif