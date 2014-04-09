#ifndef __FID_CLASS_H__
#define __FID_CLASS_H__

//--- STD Includes ----------------------------------------------------------//
#include <iostream>
//#include <fstream>
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
#include "TGraph.h"
//#include "TCanvas.h"
#include "TF1.h"

typedef vector<double> vec;

// This library consists of several frequency exaction and analysis methods
// for FIDs as well as a class to encapsulate all the ideas

namespace fid
{
	class FID {

	  public:
	  	// ctor
	  	FID(const vec& wf, const vec& tm);
	  	//FID(vec& wf, double tStart, double dt, double tZero, double tTotal);

	  	// dtor
	  	~FID();

			// frequency extraction methods
			double CalcZeroCountFreq();
			double CalcCentroidFreq();
			double CalcAnalyticalFreq();
			double CalcLorentzianFreq();
			double CalcSoftLorentzianFreq();
			double CalcExponentialFreq();
			double CalcPhaseFreq(int n=1);
			double CalcSinusoidFreq();

			// accessors
			vec& wf() {return wf_;};
			vec& tm() {return tm_;};
			const vec& power() {return power_;};
			const vec& freq() {return freq_;};
			const vec& phase() {return phase_;}
			const vec& env() {return env_;};
			const double& chi2() {return chi2_;};
			const TGraph& gr_time_series() {return gr_time_series_;};
			const TGraph& gr_freq_series() {return gr_freq_series_;};
			const TF1&    f_fit() {return f_fit_;};

		private:

			// Member Variables
			// constants
			const int kZCWidth = 40;
			const int kFitWidth = 20;
			const int kEdgeIgnore = 10;
			const double kStartThresh = 20.0;
			const double kZCAlpha = 0.8;
			const double kLowPassFreq = 100.0; // 100 kHz
			const double kMaxPhaseJump = 0.3 * 2 * M_PI;
			const double kCentroidThresh = 0.01;
			const double kTau = 2 * M_PI;
			const double kHystThresh = 0.3;

			// mutable variables
			int i_wf_; // start and stop of relevant data
			int f_wf_;
			int i_fft_;
			int f_fft_;
			double i_tm_;
			double f_tm_;
			double noise_;
			double chi2_; // Store the most recent chi2
			vec guess_;
			TF1 f_fit_;
			TGraph gr_time_series_;
			TGraph gr_freq_series_;

			// bigger data arrays
			vec wf_;
			vec tm_;
			vec power_;
			vec env_;
			vec phase_;
			vec freq_;
			vec temp_; // for random transformations

			// internal utility functions
			void CalcNoise();			
			void FindFidRange();
			void CalcPowerEnvAndPhase();
			void CalcFftFreq();
			void GuessFitParams();
			void FreqFit(TF1& func);

			// Simple n-order Butterworth filter
			inline double LowPassFilter(double f, double f0, int n){
				double denom = 1.0 + pow(f / f0, 2 * n);
				return sqrt(1.0 / denom);
			}

			// Questionable fnctions
			void GetFftPower(vec& power, vec& wf);
			void GetFftFreq(vec& freq, const vec& tm);
			void GetFftFreq(vec& freq, const int N, const double dt);
			void GetFidPhase(vec& phase, vec& wf_re);
			void GetFidEnvelope(vec& env, vec& wf_re);

	}; // FID

} // fid

#endif