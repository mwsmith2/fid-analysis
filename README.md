fid-analysis
============

A library that performs digital signal processing on free induction decay (FID) signals.  It has both cpp and python libraries.

config
======

./runtime/.default_fid_params.json is read in for default configuration 
parameters.  One should not change parameters in this.  The program also
checks for parameters in ./runtime/fid_params.json after loading the defaults.
"fid_params.h" is present in all of the project.  It holds the parameter
namespace "fid" as well as a typedef for vec=std::vector<double>.

library overview
================

### "fid_class.h"

#### class FID(vec& wf, vec& tm):
-Instantiate the class with a vector of waveform data and time data.  

Frequency extraction methods:
* .CalcZeroCountFreq()
* .CalcCentroidFreq()
* .CalcAnalyticalFreq()
* .CalcLorentzianFreq()
* .CalcSoftLorentzianFreq()
* .CalcExponentialFreq()
* .CalcPhaseFreq(int poln=1)
* .CalcSinusoidFreq()

Data accessors:
* .wf() - the waveform data
* .tm() - the time data
* .power() - the spectral density
* .freq() - frequencies associated with FFT calculations
* .phase() - the unwrapped phase
* .env() - the envelope function
* .chi2() - chi-squared from most recent fit
* .gr_time_series() - a TGraph used for FID time series (sine fit, phase fit)
* .gr_freq_series() - a TGraph used for FID spectral peak fits
* .f_fit() - TF1 function used in most recent fit


### "fid_sim.h"

#### class FidFactory()
- Instantiated with parameters specified in config file

Public Methods:
* SimulateFid(vec& wf, vec& tm)
* IdealFid(vec& wf, vec& tm)

#### class GradientFidFactory()



### "fid.h"
Meant as an alternative to create an entire FID class object.

Standalone frequency extraction methods:
* zero_count_freq(vec& wf, vec& tm)
* centroid_freq(vec& wf, vec& tm)
* analytical_freq(vec& wf, vec& tm)
* lorentzian_freq(vec& wf, vec& tm)
* soft_lorentzian_freq(vec& wf, vec& tm)
* exponential_freq(vec& wf, vec& tm)
* phase_freq(vec& wf, vec& tm, int n=1)
* sinusoid_freq(vec& wf, vec& tm)

FFT helper functions:
* fft_power(vec& power, vec& wf)
* fft_freq(vec& freq, vec& tm)
* fft_freq(vec& freq, int N, double dt)
* fid_phase(vec& phase, vec& wf_re)
* fid_envelope(vec& env, vec& wf_re)
* ideal_fid(vec& wf, vec& tm, double f, double phi=0.0, double snr=100.0, 			double tau=10.0, double t0=0.0)

### "fid_utilities.h"
A bunch of functions that may come in handy.

* DrawFID(vec &wf, vec &tm, string filename, string title);
* DrawFID(FID &my_fid, string filename, string title);

* AddWhiteNoise(vec &wf, double snr=100.0);
  - Add Gaussian noise to the given waveform

* ConstructTimeVector(int num_times, double t0, double dt, vec &tm);
  - Get a time vector for the FID

* ConstructLinearGradient(int num_points, vec &grad);
  - Get a linear gradient with a max(abs(grad)) == 1 and mean(grad) == 0

* ConstructQuadraticGradient(int num_points, vec &grad);
  - Get a quadratic gradient with a max(abs(grad)) == 1 and mean(grad) = 0

* ConstructSweepRange(vector<T>& range_vec){
  - Use for getting sweep ranges from parameter namespace
