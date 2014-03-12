#ifndef __FID_H__
#define __FID_H__

//--- STD Includes ----------------------------------------------------------//
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

//--- Other Includes --------------------------------------------------------//
#include <fftw3.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"

typedef vector<double> vec;

// This library consists of several frequency exaction and analysis methods
// for FIDs as well as a class to encapsulate all the ideas

namespace fid
{
	// Define some constants
	const int fit_w = 20;
	const double start_thresh = 20;
	const int zc_width = 40;
	const double zc_alpha = 0.5;
	const double ph_max_jump = 0.7;
	const int ph_edge_ignore = 10;

	// Declare all the frequency extraction methods as standalone functions
	double calcZeroCountFreq(vec& wf, const vec& tm);
	double calcCentroidFreq(vec& wf, const vec& tm);
	double calcAnalyticalFreq(vec& wf, const vec& tm);
	double calcLorentzianFreq(vec& wf, const vec& tm);
	double calcSoftLorentzianFreq(vec& wf, const vec& tm);
	double calcExponentialFreq(vec& wf, const vec& tm);
	double calcPhaseFreq(vec& wf, const vec& tm, int n=1);
	double calcSinusoidFreq(vec& wf, const vec& tm);

	// utility functions
	void getFFTPower(vec& power, vec& wf);
	void getFFTFreq(vec& freq, const vec& tm);
	void getFFTFreq(vec& freq, const int N, const double dt);
	void getFIDPhase(vec& phase, vec& wf_re);
	void getFIDEnvelope(vec& env, vec& wf_re);
	void getIdealFID(vec& wf, vec& tm, double f=12.345678, double phi=0.0, 
		double s2n=40.0, double tau=10.0, double t0=0.0);

} // fid

#endif