#ifndef FID_INCLUDE_FID_H_
#define FID_INCLUDE_FID_H_

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
#include "TGraph.h"
#include "TF1.h"

//--- project includes ------------------------------------------------------//
#include "fid_params.h"
#include "fid_math.h"

// This library consists of several frequency exaction and analysis methods
// for FIDs as well as a class to encapsulate all the ideas

namespace fid
{
	// Declare all the frequency extraction methods as standalone functions
	double zero_count_freq(vec& wf, const vec& tm);
	double spd_centroid_freq(vec& wf, const vec& tm);
	double analytical_freq(vec& wf, const vec& tm);
	double lorentzian_freq(vec& wf, const vec& tm);
	double soft_lorentzian_freq(vec& wf, const vec& tm);
	double exponential_freq(vec& wf, const vec& tm);
	double phase_freq(vec& wf, const vec& tm, int n=1);
	double sinusoid_freq(vec& wf, const vec& tm);

	void ideal_fid(vec& wf, vec& tm, double f, double phi=0.0, double snr=100.0, double tau=10.0, double t0=0.0);

} // fid

#endif