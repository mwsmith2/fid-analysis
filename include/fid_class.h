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
			const vec& power() {return power_;};
			const vec& freq() {return freq_;};
			const vec& phase() {return phase_;}
			const vec& env() {return env_;};

		private:

			// Member Variables
			// constants
			const int kZCWidth = 40;
			const int kFitWidth = 20;
			const int kEdgeIgnore = 10;
			const double kStartThresh = 20.0;
			const double kZCAlpha = 0.5;
			const double kMaxPhaseJump = 0.7 * 2 * M_PI;
			const double kCentroidThresh = 0.1;
			const double kTau = 2 * M_PI;

			// mutable variables
			int i_wf_; // start and stop of relevant data
			int f_wf_;
			int i_fft_;
			int f_fft_;
			double i_tm_;
			double f_tm_;
			double noise_;
			vec guess_;
			TGraph gr_peak_;

			// bigger data arrays
			vec wf_;
			vec tm_;
			vec power_;
			vec env_;
			vec phase_;
			vec freq_;
			vec temp_; // for random transformations

			// utility functions
			void CalcNoise();			
			void FindFidRange();
			void CalcPowerEnvAndPhase();
			void CalcFftFreq();
			void GuessFitParams();
			void FreqFit(TF1& func);

			// Questionable fnctions
			void GetFftPower(vec& power, vec& wf);
			void GetFftFreq(vec& freq, const vec& tm);
			void GetFftFreq(vec& freq, const int N, const double dt);
			void GetFidPhase(vec& phase, vec& wf_re);
			void GetFidEnvelope(vec& env, vec& wf_re);

	}; // FID

} // fid

#endif