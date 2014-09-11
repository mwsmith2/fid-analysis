#include "fid.h"

namespace fid
{
	// Find the frequency by counting zeros and interpolating time 
	// at the first and last zeros.
	double zero_count_freq(vec& wf, const vec& tm)
	{
		// Find the noise level at the head and tail.
		double head = stdev(wf_.begin(), wf_.begin() + params::zc_width);
		double tail = stdev(wf_.rbegin(), wf_.rbegin() + params::zc_width);

		// Take the smaller of the two.
		double noise = (tail < head) ? (tail / w) : (head / w);

		// And set the rms noise.
		noise = std::sqrt(noise); // Want the RMS

		// Now find the starting and ending points
		double thresh = 20 * noise;
		auto it_i = std::find_if(wf.begin(), wf.end(), 
				[thresh](double x){return std::abs(x) > thresh;});

		auto it_f = std::find_if(wf.rbegin(), wf.rend(), 
				[thresh](double x){return std::abs(x) > thresh;});

		int fid_i = std::distance(wf.begin(), it_i);
		int fid_f = std::distance(it_f, wf.rend());

		// set up vectors to hold relevant stuff about the important part
		vec wf_smooth (fid_f - fid_i, 0.0);
		vector<int> sign_vec (fid_f - fid_i, 0);
		vector<int> diff_vec (fid_f - fid_i, 0);

		// smooth the waveform, exponential moving average algorithm
		double a = params::zc_alpha;
		std::partial_sum(wf.begin() + fid_i, wf.begin() + fid_f, 
			wf_smooth.begin(), 
			[a](double sum, double x) {return sum * a + x * (1.0 - a);});

		// Get a vector of just the signs
		std::transform(wf_smooth.begin(), wf_smooth.end(), sign_vec.begin(),
			[](double x) {return x >= 0.0 ? 1 : -1;});

		// Find out when the sign changes
		std::adjacent_difference(sign_vec.begin(), sign_vec.end(), 
			diff_vec.begin());

		// Count the sign changes
		int zeros = std::count_if(diff_vec.begin()+1, diff_vec.end(), 
			[](int x) {return x != 0;});

		// Use linear interpolation to find the zeros more accurately
		auto it_i2 = std::find_if(diff_vec.begin()+1, diff_vec.end(),
			[](int x) {return x != 0;});
		auto it_f2 = std::find_if(diff_vec.rbegin(), diff_vec.rend(),
			[](int x) {return x != 0;});

		int i = fid_i + std::distance(diff_vec.begin(), it_i2) - 1;
		int f = fid_i + std::distance(it_f2, diff_vec.rend()) - 1;

		double frac = std::abs(wf[i] / (wf[i-1] - wf[i]));
		double ti = frac * tm[i-1] + (1.0 - frac) * tm[i];

		frac = std::abs(wf[f] / (wf[f-1] - wf[f]));
		double tf = frac * tm[f-1] + (1.0 - frac) * tm[f];

		return 0.5 * (zeros - 1) / (tf - ti);
	}

	double centroid_freq(vec& wf, const vec& tm)
	{
		// Get the Spectral Density and frequencies
		vec power;
		vec freq;
		fft_power(power, wf);
		fft_freq(freq, tm);

		// Find the peak power
		double max = *std::max_element(power.begin(), power.end());

		// Find the half-max indices
		int it_i = std::distance(power.begin(), 
			std::find_if(power.begin(), power.end(), 
				[max](double x) {return x > 0.1 * max;}));

		// reverse the iterators
		int it_f = -1 * std::distance(power.rend(),
			std::find_if(power.rbegin(), power.rend(), 
				[max](double x) {return x > 0.1 * max;}));

		// Now compute the power weighted average
		double pwfreq = 0.0;
		double pwsum = 0.0;

		for (int i = it_i; i < it_f; i++){
			pwfreq += power[i] * freq[i];
			pwsum  += power[i];
		}

		return pwfreq / pwsum;
	}

	double analytical_freq(vec& wf, const vec& tm)
	{
		// @todo check the algebra on this guy
		// Set the fit function
		std::string num1("[2] * ([0]^2 - 0.5 * [1] * [0] * sin(2 *[4])");
		std::string num2(" + (0.5 * [1])^2 + x^2 - [0]^2 * sin([4])^2) / ");
		std::string den1("((0.5 * [1])^2 + 2 * (x^2 - [0]^2) * (x^2 + [0]^2)");
		std::string den2(" + (x^2 - [0]^2)^2) + [3]");

		TF1 f_fit = TF1("f_fit", (num1 + num2 + den1 + den2).c_str());

		// @bug these should be calculated not, hardcoded
		// Guess the parameters
		f_fit.SetParameter(0, 23);
		f_fit.SetParameter(1, 0.05);
		f_fit.SetParameter(2, 350);
		f_fit.SetParameter(3, 0.0);
		f_fit.SetParameter(4, 0.0);
		f_fit.SetParLimits(4, 0.0, 6.29);

		// Get the Spectral Density and frequencies
		vec power;
		vec freq;
		fft_power(power, wf);
		fft_freq(freq, tm);

		// Make a TGraph to fit
		TGraph gr_fit = TGraph(freq.size(), &freq[0], &power[0]);

		int max_idx = std::distance(power.begin(), 
			std::max_element(power.begin(), power.end()));

		int w = params::fit_width;
		gr_fit.Fit("f_fit", "QRNM", "", freq[max_idx - w], freq[max_idx + w]);

		return f_fit.GetParameter(0);
	}

	double lorentzian_freq(vec& wf, const vec& tm)
	{
		// Set the fit function
		TF1 f_fit = TF1("f_fit", 
			"[2] / (1 + ((x - [0]) / (0.5 * [1]))^2) + [3]");

		// @bug these should be calculated not, hardcoded
		// Guess the parameters
		f_fit.SetParameter(0, 23);
		f_fit.SetParameter(1, 0.05);
		f_fit.SetParameter(2, 350);
		f_fit.SetParameter(3, 0.0);
		f_fit.SetParameter(4, 2.0);
		f_fit.SetParLimits(4, 1.0, 3.0);

		// Get the Spectral Density and frequencies
		vec power;
		vec freq;
		fft_power(power, wf);
		fft_freq(freq, tm);

		// Make a TGraph to fit
		TGraph gr_fit = TGraph(freq.size(), &freq[0], &power[0]);

		int max_idx = std::distance(power.begin(), 
			std::max_element(power.begin(), power.end()));

		int w = params::fit_width;
		gr_fit.Fit("f_fit", "QRNM", "", freq[max_idx - w], freq[max_idx + w]);

		return f_fit.GetParameter(0);
	}

	double soft_lorentzian_freq(vec& wf, const vec& tm)
	{
		// Set the fit function
		TF1 f_fit = TF1("f_fit", 
			"[2] / (1 + ((x - [0]) / (0.5 * [1]))^[4]) + [3]");

		// @bug these should be calculated not, hardcoded
		// Guess the parameters
		f_fit.SetParameter(0, 23);
		f_fit.SetParameter(1, 0.05);
		f_fit.SetParameter(2, 350);
		f_fit.SetParameter(3, 0.0);
		f_fit.SetParameter(4, 2.0);
		f_fit.SetParLimits(4, 1.0, 3.0);

		// Get the Spectral Density and frequencies
		vec power;
		vec freq;
		fft_power(power, wf);
		fft_freq(freq, tm);

		// Make a TGraph to fit
		TGraph gr_fit = TGraph(freq.size(), &freq[0], &power[0]);

		int max_idx = std::distance(power.begin(), 
			std::max_element(power.begin(), power.end()));

		int w = params::fit_width;
		gr_fit.Fit("f_fit", "QRNM", "", freq[max_idx - w], freq[max_idx + w]);

		return f_fit.GetParameter(0);
	}

	double exponential_freq(vec& wf, const vec& tm)
	{
		// Set the fit function
		TF1 f_fit = TF1("f_fit", "[2] * exp(-abs(x - [0]) / [1]) + [3]");

		// Get the Spectral Density and frequencies
		vec power;
		vec freq;
		fft_power(power, wf);
		fft_freq(freq, tm);

		// Make a TGraph to fit
		TGraph gr_fit = TGraph(freq.size(), &freq[0], &power[0]);

		// @bug these should be calculated not, hardcoded
		// Guess the parameters
		f_fit.SetParameter(0, 23);
		f_fit.SetParameter(1, 0.05);
		f_fit.SetParameter(2, 350);
		f_fit.SetParameter(3, 0.0);

		int max_idx = std::distance(power.begin(), 
			std::max_element(power.begin(), power.end()));

		int w = params::fit_width;
		gr_fit.Fit("f_fit", "QRNM", "", freq[max_idx - w], freq[max_idx + w]);

		return f_fit.GetParameter(0);	
	}

	double phase_freq(vec& wf, const vec& tm, int n)
	{
		// Container for the phase
		vec phase;

		// Get the phase
		fid_phase(phase, wf);

		// Make a TGraph to fit
		TGraph gr_fit = TGraph(tm.size(), &tm[0], &phase[0]);

		// Now set up the polynomial phase fit
		char fcn[20];
		sprintf(fcn, "pol%d", n);
		TF1 f_fit = TF1("f_fit", fcn);

		// Guess the parameters
		f_fit.SetParameter(1, 23 * 2 * M_PI);

		// Find the fit range
		// Need the noise level
		double w = params::fit_width;
		double noise = std::accumulate(wf.begin(), wf.begin() + w, 0.0, 
			[](double x, double y){return x + y * y;}); 
		double temp = std::accumulate(wf.rbegin(), wf.rbegin() + w, 0.0, 
			[](double x, double y){return x + y * y;}); 

		noise = (temp < noise) ? (temp / w) : (noise / w);
		noise = std::pow(noise, 0.5); // Take the 

		// Now find the starting and ending points
		double thresh = params::start_thresh * noise;
		auto it_i = std::find_if(wf.begin(), wf.end(), 
				[thresh](double x){return std::abs(x) > thresh;});

		auto it_f = std::find_if(wf.rbegin(), wf.rend(), 
				[thresh](double x){return std::abs(x) > thresh;});

		int fid_i = std::distance(wf.begin(), it_i);
		int fid_f = std::distance(it_f, wf.rend());

		// Adjust to ignore the edges
		fid_i += params::edge_ignore;
		fid_f -= params::edge_ignore;

		// Do the fit.
		gr_fit.Fit("f_fit", "QRNM", "", tm[fid_i], tm[fid_f]);

		return 0.5 * f_fit.GetParameter(1) / M_PI;
	}

	double sinusoid_freq(vec& wf, const vec& tm)
	{
		// Get the FID envelope function
		vec env;
		fid_envelope(env, wf);

		// Normalize the waveform by the envelope
		vec wf_nm(wf.size(), 0.0);

		std::transform(wf.begin(), wf.end(), env.begin(), wf_nm.begin(),
			[](double x_wf, double wf_env) {return x_wf / wf_env;} );

		TGraph gr_fit = TGraph(wf.size(), &tm[0], &wf_nm[0]);

		TF1 f_fit = TF1("f_fit", "[1] * sin([0] * x + [2])");
		f_fit.SetParameters(23.0 * 2 * M_PI, 1.0, 0.0);


		// Find the fit range, Need the noise level.
		double w = params::fit_width;
		double noise = std::accumulate(wf.begin(), wf.begin() + w, 0.0, 
			[](double x, double y){return x + y * y;}); 
		double temp = std::accumulate(wf.rbegin(), wf.rbegin() + w, 0.0, 
			[](double x, double y){return x + y * y;}); 

		noise = (temp < noise) ? (temp / w) : (noise / w);
		noise = std::pow(noise, 0.5); // Take the 

		// Now find the starting and ending points
		double thresh = params::start_thresh * noise;
		auto it_i = std::find_if(wf.begin(), wf.end(), 
				[thresh](double x){return std::abs(x) > thresh;});

		auto it_f = std::find_if(wf.rbegin(), wf.rend(), 
				[thresh](double x){return std::abs(x) > thresh;});

		int fid_i = std::distance(wf.begin(), it_i);
		int fid_f = std::distance(it_f, wf.rend());

		// Adjust to ignore the edges
		fid_i += params::edge_ignore;
		fid_f -= params::edge_ignore;

		// Do the fit.
		gr_fit.Fit("f_fit", "QRNM", "", tm[fid_i], tm[fid_f]);

		return 0.5 * f_fit.GetParameter(0) / M_PI;
	}

	// FFT Utility Functions

	// Wrapper for a simple 1D fft from fftw
	void fft_power(vec& power, vec& wf){

		// Set up the FFT plan
		int N = wf.size();
		double Nroot = std::sqrt(N);

		fftw_complex *fft;
		fft = new fftw_complex[N];

		fftw_plan wf_to_fft;
		wf_to_fft = fftw_plan_dft_r2c_1d(N, &wf[0], &fft[0], FFTW_ESTIMATE);
		fftw_execute(wf_to_fft);

		power.reserve(N);
		power.resize(0);

		int n = 0;
		if (N % 2 == 0){
			n = N / 2 + 1;
		} else {
			n = (N + 1) / 2;
		}

		for (int i = 0; i < n; i++){

			static double rp;
			static double ip;
			rp = fft[i][0];
			ip = fft[i][1];

			power.push_back((rp*rp + ip*ip) / Nroot);
		}

		fftw_destroy_plan(wf_to_fft);
		delete[] fft;
		return;
	}

	// Helper function to get frequencies for FFT
	void fft_freq(vec& freq, const vec& tm){

		int N = tm.size();
		double df = (N - 1) / (tm[N-1] - tm[0]);

		if (tm.size() % 2 == 0){

			freq.reserve(N / 2 + 1);
			freq.resize(0);
			
			for (int i = 0; i < N / 2 + 1; i++){
				freq.push_back(i * df / N);
			}

		} else {

			freq.reserve((N + 1) / 2);
			freq.resize(0);

			for (int i = 0; i < (N + 1) / 2 + 1; i++){
				freq.push_back(i * df / N);
			}
		}

		return;
	}

	void fid_phase(vec& ph, vec& wf_re){

		// Set up the FFT plan
		int N = wf_re.size();
		int n = N/2 + 1;
		double Nroot = std::sqrt(N);

		fftw_complex *fft;
		fft = new fftw_complex[n]; 

		fftw_plan wf_to_fft;
		wf_to_fft = fftw_plan_dft_r2c_1d(N, &wf_re[0], &fft[0], FFTW_ESTIMATE);
		fftw_execute(wf_to_fft);

		// Now get ready for Hilbert transform
		vec wf_im(wf_re.size(), 0.0);
		double temp;

		// Flip the phase
		for (int i = 0; i < n; i++){
			temp = fft[i][0];
			fft[i][0] = fft[i][1];
			fft[i][1] = -1 * temp;
		}

		fftw_plan fft_to_wf;
		fft_to_wf = fftw_plan_dft_c2r_1d(N, &fft[0], &wf_im[0], FFTW_ESTIMATE);
		fftw_execute(fft_to_wf);

		// fftw is unnormalized fft, so rescale
		for (auto it = wf_im.begin(); it != wf_im.end(); it++){
			*it /= N;
		}

		// Now compute the phase
		ph.resize(wf_re.size());
		std::transform(wf_re.begin(), wf_re.end(), wf_im.begin(), ph.begin(),
			[](double re, double im) {return std::atan2(im, re);});
	
		// Need to find the point to start unwrapping the phase.
		// Determine the noise first
		double w = params::fit_width;
		double noise = std::accumulate(wf_re.begin(), wf_re.begin() + w, 0.0, 
			[](double x, double y){return x + y * y;}); 
		temp = std::accumulate(wf_re.rbegin(), wf_re.rbegin() + w, 0.0, 
			[](double x, double y){return x + y * y;}); 

		noise = (temp < noise) ? (temp / w) : (noise / w);
		noise = std::pow(noise, 0.5); // Take the 

		// Now find the starting and ending points
		double thresh = params::start_thresh * noise;
		auto it_i = std::find_if(wf_re.begin(), wf_re.end(), 
				[thresh](double x){return std::abs(x) > thresh;});

		int fid_i = std::distance(wf_re.begin(), it_i);

		// Now unwrap the phase
		for (auto it = ph.begin() + fid_i + 1; it != ph.end(); it++){
			static int k;
			if (it == ph.begin()) k = 0; // reset

			*it += 2 * k * M_PI;

			// Check for jumps, both positive and negative.
			if (*(it-1) - *it > 2 * params::max_phase_jump * M_PI) {

				k++;
				*it += 2 * M_PI;

			} else if (*(it-1) - *it < -2 * params::max_phase_jump * M_PI){

				k--;
				*it -= 2 * M_PI;

			}
		}

		delete[] fft;
		fftw_destroy_plan(wf_to_fft);
		fftw_destroy_plan(fft_to_wf);

		return;
	}

	void fid_envelope(vec& env, vec& wf_re)
	{

		// Set up the FFT plan
		int N = wf_re.size(); // size of data
		int n = N/2 + 1;      // size of rfft
		double Nroot = std::sqrt(N);

		fftw_complex *fft = new fftw_complex[n]; 

		fftw_plan wf_to_fft;
		wf_to_fft = fftw_plan_dft_r2c_1d(N, &wf_re[0], &fft[0], FFTW_ESTIMATE);
		fftw_execute(wf_to_fft);

		// Now get ready for Hilbert transform
		vec wf_im(wf_re.size(), 0.0);
		double temp;

		// Flip the phase
		for (int i = 0; i < n; i++){
			temp = fft[i][0];
			fft[i][0] = fft[i][1]; 
			fft[i][1] = -1 * temp;
		}

		fftw_plan fft_to_wf;
		fft_to_wf = fftw_plan_dft_c2r_1d(N, &fft[0], &wf_im[0], FFTW_ESTIMATE);
		fftw_execute(fft_to_wf);

		// fftw is unnormalized fft, so rescale
		for (auto it = wf_im.begin(); it != wf_im.end(); it++){
			*it /= N;
		}

		// Now fill the envelope function
		env.resize(wf_re.size());
		std::transform(wf_re.begin(), wf_re.end(), wf_im.begin(), env.begin(),
			[](double re, double im) {return std::sqrt(im * im + re * re);});

		// Clean up
		delete[] fft;
		fftw_destroy_plan(wf_to_fft);
		fftw_destroy_plan(fft_to_wf);

		return;
	}

	// Generate an ideal FID
	void ideal_fid(vec& wf, vec& tm, double f, double phi, 
		double snr, double tau, double t0){

		wf.reserve(tm.size());
		wf.resize(0);

		// Define the waveform
		double temp;

		for (auto it = tm.begin(); it != tm.end(); it++){

			if (*it >= t0){
				temp = std::exp(-(*it - t0) / tau);
				temp *= std::sin((*it) * 2 * M_PI * f + phi);
				wf.push_back(temp);

			} else {
				wf.push_back(0.0);

			}
		} 

		// Add some noise
		static std::default_random_engine gen(0);
		std::normal_distribution<double> norm(0.0, 1.0 / snr);
		for (auto it = wf.begin(); it != wf.end(); it++){
			*it += norm(gen);
		}

		return;
	}

}



