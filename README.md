# fid-analysis

A library that performs digital signal processing on free induction decay (FID) signals.  It has a cpp library and the intention of adding a python wrapper for the library.

## Getting Started

### Platforms
The library builds on both Linux and OS X.  I don't think there is any reason it can't build on Windows, but the Makefile is not configured to facilitate building on Windows.

### Prerequsites
The fid-analysis library depends on [fftw3](http://www.fftw.org), [armadillo](http://arma.sourceforge.net), [boost](https://boost.org), [ROOT](https://root.cern.ch/), and several c++11 features.

#### Linux

See the [ROOT website](https://root.cern.ch/downloading-root) for instructions on getting ROOT set up.

On Ubuntu
```bash
apt-get install fftw-dev armadillo-dev boost-dev
```

On Fedora (or CentOS, Scientific Linux)
```bash
yum install fftw-devel armadillo-devel boost-devel
```

#### OS X

I highly recommend setting up [homebrew](http://brew.sh/) if you aren't using it.  You can even manage ROOT with the science tap.

```bash
brew install fftw boost armadillo root
```

### Installation
First, you'll want to download the repository.  Then, proceed to build and install it.

```
git clone https://github.com/mwsmith2/fid-analysis
cd fid-analysis/cpp
make install
```

To use the library you just need to include fid.h and use the linker command -lfid when compiling.

### Quick Example
Example files can be found in the examples folder of the repository along with a Makefile to compile them.  Let's look at one right here.

```cpp
#include <iostream>
#include "fid.h"

using namespace fid;

int main(int argc, char **argv)
{
  FID myfid(argv[1]);

  std::cout << "Frequency, Error" << std::endl;
  std::cout << myfid.GetFreq() << ", " << myfid.GetFreqError() << std::endl;

  return 0;
}
```

The file in this case is a text file with the format:

```
<time 0> <amp 0>
<time 1> <amp 1>
...
```

## Library Overview
The library is fairly small, consisting of a few classes and some utility functions.

### FID Class

#### Constructors

Upon instantiation the FID class performs cursory analysis and computes the fft, spectral density, phase, and envelope. There are several ways to pass the data into the class.

* FID(std::string fidfile)
	- read FID data in from a text file
	- format is `<time 0> <amp 0>`

* FID(std::vector<double> wf, std::vector<double> tm)
	- the input vectors data is set beforehand

* FID(std::vector<double> wf)
	- the input times are set to integer values for each datum in wf

#### Frequency Extraction Methods

All frequency extraction techniques can be called by setting the current class method and calling GetFreq(<method-string>) or by an explicit function call.

* Zero Counting
	- GetFreq("ZC")
	- GetFreq("ZEROCOUNT")
	- CalcZeroCountFreq()
* Spectral Centroid
	- GetFreq("CN")
	- GetFreq("CENTROID")
	- CalcCentroidFreq()
* Lorentzian Peak Fit
	- GetFreq("LZ")
	- GetFreq("LORENTZIAN")
	- CalcLorentzianFreq()
* Exponential Peak Fit
	- GetFreq("EX")
	- GetFreq("EXPONENTIAL")
	- CalcExponentialFreq()
* Analytical Peak Fit
	- GetFreq("AN")
	- GetFreq("ANALYTICAL")
	- CalcAnalyticalFreq()
* Polynomial Phase Fit
	- GetFreq("PH")
	- GetFreq("PHASE")
	- CalcPhaseFreq(int poln=1)
* Normalized Sine Fit
	- GetFreq("SN")
	- GetFreq("SINUSOID")
	- CalcSinusoidFreq()

#### Accessors
Accessors return a const reference to several internal class variables.

* wf() - waveform vector
* tm() - time vector
* power() - pectral density
* freq() - frequencies associated with FFT calculation
* phase() - unwrapped phase from hilbert transform
* env() - envelope function
* chi2() - chi-squared from most recent fit
* f_fit() - TF1 function used in most recent fit
* gr_time_series() - a TGraph used for FID time series (sine fit, phase fit)
* gr_freq_series() - a TGraph used for FID spectral peak fits

### FastFid Class

This is a trimmed version of the FID class for use when computation time matters, for instance online analysis in data acqusition.  It doesn't do any frequency analysis calculations (no FFT), and only does zero crossing.

### FidFactory Class

An FID simulation class that uses the idealized FID form and Bloch equations to create FID waveform.  The types of FIDs are detailed below.  Each method has option to add noise or not

#### Sim Params

The way to set parameters for simulated FIDs is by setting certain parameters in the fid::sim namespace.

```c++
int seed = 0; // seed for random noise generator
double dt_integration = 2.0e-5; // integration step size
double snr = 90000.0; // signal to noise ratio for adding noise
double amplitude = 2000.0; // half peak to peak amplitude for FID
double baseline = 21000.0; // baseline offset
int num_samples = 10000; // numbers of samples in final FID
double start_time = -1.0; // t0 for FID
double delta_time = 0.001; // time step
double freq_ref = 950.0; // mixdown clock frequency
double freq_larmor = 997.0; // nmr frequency for 
double freq_cut_ratio = 0.1; // sets frequency for lowpass filter
double mixdown_phi = 950.0; // starting phase for mixdown clock
std::vector<double> spin_0 = {0.0, 0.0, 1.0}; // intial spin for sims
double gamma_1 = 0.05; // spin relaxation time
double gamma_2 = 0.05; // decoherence relaxation time
double gamma_g = 1.0;  // proton gyromagnetic ratio
double omega_r = 50.0; // RF frequency for NMR pulser
double t_pulse = 0.005; // time of NMR pulser kick
```

#### Ideal FIDs

The simplest model is just an exponential decay multiplied with a sine function.

```c++
void IdealFid(std::vector<double>& wf, 
	          std::vector<double>& tm, 
		      bool withnoise=false,
		 	  bool discretize=false)
```

#### Bloch Equation FIDs

Integrate an FID based on parameters from fid::sim.

```c++
void SimulateFid(std::vector<double>& wf, 
                 std::vector<double>& tm, 
                 bool withnoise=false,
                 bool discretize=false);
```


#### Gradient FIDs

Does a summation of several FIDs to simulate the effect that a small gradient would have along the active region of an NMR probe.  The caveat here is that you need to generate a ROOT file filled with previous simulated FIDs (I may force it to simulate FIDs with the Bloch Equations if not found but not yet implemented).  The default place it looks is `~/.fid/sim_fids.root`, but this can be changed using the variables `fid::grad::root_file`.

```c++
GradientFid(std::vector<double>& gradient, 
            std::vector<double>& wf, 
            bool withnoise=false,
            bool discretize=false);
```

### FID Math Functions

Several of the core math functions are implemented outside of classes, so that they can be used in a general way. Most of them are related to FFT and frequency analysis.

```c++
std::vector<std::complex<double>> dsp::fft(const std::vector<double> &wf)

std::vector<double> dsp::ifft(const std::vector<std::complex<double>>& fft)

std::vector<double> dsp::psd(const std::vector<double>& wf)

std::vector<double> dsp::fftfreq(const std::vector<double>& tm)

std::vector<double> dsp::hilbert(const std::vector<double>& wf)

std::vector<double> dsp::phase(const std::vector<double>& wf)

std::vector<double> dsp::envelope(const std::vector<double>& wf)

std::vector<double> dsp::savgol3(const std::vector<double>& wf)

std::vector<double> dsp::savgol5(const std::vector<double>& wf)
```

The fft and inverse fft (ifft) make the assumption that the original waveform was purely real.  The Hilbert transform, phase, and envelope functions are computed by assuming that the signal is harmonic.  The last two functions (`savgol3, savgol5`) are Sovitsky-Golay filters which compute the derivative of data in a way that is robust against noisy data.

### Utility Functions

#### Configuration

The library uses several internal parameters from the namespace fid which can be set using an external JSON file.  The function to call is just `fid::load_params(string filename)`.  Each variable is set by using it's internal name/namespace structure.  An example config file is given in `cpp/config/user_fid_params.json`.  If variables are not present in the config file, they retain their default values defined in `cpp/include/fid_params_extdef.h`.


#### Saving Fit Plots

The functions that let one save fit results and residuals (very useful for debugging) are the following:

```c++
void draw_fid_time_fit(const FID &my_fid, 
                       std::string fname, 
                       std::string title);

void draw_fid_freq_fit(const FID &my_fid, 
                       std::string fname, 
                       std::string title);

void draw_fid_time_res(const FID &my_fid, 
                       std::string fname, 
                       std::string title);

void draw_fid_freq_res(const FID &my_fid, 
                       std::string fname, 
                       std::string title);
```

#### Saving FID Data

There are a few functions along these lines.  One saves the FID data itself to a text file.  Another saves frequency extraction results to a file stream.  Yet another save frequency extraction results for different types of phase fits to a file stream.

```c++
void write_fid_file(std::string fname, 
                    const std::vector<double> &wf, 
                    const std::vector<double> &tm);

void calc_freq_write_csv(FID &my_fid, std::ofstream &out);

void calc_phase_freq_write_csv(FID &my_fid, std::ofstream &out);
```

## Contributing

Feature requests can be made to the github page, and if you would like to contribute just email me at mwsmith2@uw.edu.