"""Analyzer class for FID waveforms."""

import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import splrep, sproot
from pandas import read_csv

# Define some functions will be used to fit peaks
def real_spectral_density(x, x0, gamma, a, b):

	c1 = (0.5 * gamma)**2
	c2 = (x**2 - x0**2)
	c3 = (x**2 + x0**2)

	num = a * (c1 + 0.5 * (c2 + c3))
	den = c1**2 + 2 * c1 * c3 + c2**2
	return num / den + b

def lorentzian(x, x0, gamma, a, b):

	num = 0.5 * gamma * a
	den = np.pi * ((x - x0)**2 + (0.5 * gamma)**2)
	return num / den + b

def soft_lorentzian(x, x0, gamma, a, b, alpha):

	num = 0.5 * gamma * a
	den = np.pi * (abs(x - x0)**alpha + abs(0.5 * gamma)**alpha)
	return num / den + b

def exp_peak(x, x0, tau, c, b):

	return np.exp(-abs((x - x0)/tau) + c) + b

""" The FID class
It has several frequency extraction methods.

"""

class FID:
	
	def __init__(self, fid):

		# Get the data from file if its a file, else it's just two arrays.
		if not (isinstance(fid, np.ndarray)):  

			data = np.ndfromtxt(fid)
			#data = read_csv(fid_file, sep=' ', header=0).values

		else:

			data = fid

		self.tm = data[:, 0]
		self.wf = data[:, 1] - data[:, 1].mean() # center waveform
		del data

		# Compute time spacing.
		self.N  = self.tm.shape[0]
		self.dt = (self.tm[-1] - self.tm[0]) / (self.N - 1)

		self.freq_guess = 0
		self.guess = []
		self.w0 = 50

		# Find the beginning and end indexes
		self.fid_i = -1
		self.fid_f = -1

		# Estimate the noise.
		noise = self.wf[:100].std()	
		if (noise > self.wf[-100:].std()):
			noise = self.wf[-100:].std()
		self.noise = noise
		
		# Compute the fft.
		self.fft  = np.fft.rfft(self.wf, n=self.N) 
		self.fft /= self.N**0.5 # Smaller numbers for better numerics
		self.freq = np.fft.fftfreq(self.wf.shape[0], d=self.dt)[:self.N/2+1]
		self.freq[-1] *= -1 # Stupid trick since no rfftfreq helper function
		self.power = np.abs(self.fft)**2

		# Compute the envelope function
		self.wf_im = np.fft.irfft(self.fft * self.N**0.5 * -1j, n=self.N)
		self.env = np.sqrt(self.wf**2 + self.wf_im**2)

		self.phase = np.empty([0]) # Fill later

		# BUG this was tested on noiseless waveforms
		thresh = 0.15 * self.env.max()
		if (thresh > 20 * self.noise):
			thresh = 20 * self.noise

		z = np.where(np.diff(np.abs(self.wf) > thresh) != 0)[0]

		if (z.shape[0] != 0):
			self.fid_i = z[0]
		else:
			self.fid_i = 10

		if (z.shape[0] < 2):
			self.fid_f = self.N - 10 # Arbitrary
		else:
			self.fid_f = z[-1]

		# Need to interpolate to find start and end times.
		d = self.w0
		j = self.fid_i
		self.ti = sproot(splrep(self.tm[j:j+d], self.wf[j:j+d]))[0]
		j = self.fid_f
		self.tf = sproot(splrep(self.tm[j-d:j], self.wf[j-d:j]))[0]

		# Define a dictionary to access different frequency compuations
		self.method = {}
		for i in range(20):
			self.method[i] = i

		self.method['zc'] = 0 # Zero Counting
		self.method['cn'] = 1 # Centroid
		self.method['lo'] = 2 # Lorentzian Peak
		self.method['sl'] = 3 # Soft Lorentzian
		self.method['ex'] = 4 # Exponential Peak
		self.method['p1'] = 5 # 1st Order Phase Polynomial
		self.method['p2'] = 6 # 2nd Order Phase Polynomial
		self.method['p3'] = 7 # 3rd Order Phase Polynomial
		self.method['os'] = 8 # Sinusoidal with Envelop Normalization



	def getZeroCountFreq(self):

		ncross = 0
		pos = True
		thresh = 0.7 #5 * self.noise
		min_sep = 8

		# Slice out the relevant part of the FID.
		d = self.wf[self.fid_i:self.fid_f]
		e = self.env[self.fid_i:self.fid_f]

		# Now find the times where the transition through zero happens
		idx = np.where(np.abs(d) > thresh * e)[0]
		ncross = np.where(np.diff(idx) > min_sep)[0].shape[0]

		return 0.5 * (ncross) / (self.tf - self.ti)

	# @bug not implemented
	def getInfPointFreq(self):

		print "Not Implemented yet."
		return None

		ninf = 0
		i = 0
		f = 0

		der2_window = np.array([1, -2, 1])
		der2_window /= self.dt**2
		
		der2 = np.convolve(self.wf, der2_window, mode='same')
		pos = der2[self.fid_i] >= 0

		for n in range(self.fid_i, self.fid_f):

			if (pos != der2[n] >= 0):
				if (i == 0):
					i = n
				pos = der2[n] >= 0
				ninf += 1
				f = n

		freq_infl = 0.5 * ninf / (self.tm[f] - self.tm[i])

		self.der2 = der2
		return freq_infl

	def getCentroidFreq(self):

		# Compute the frequency using the FFT centroid
		nmax = self.power.argmax()
		if (nmax > self.w0):
			w = self.w0
		else:
			w = nmax

		i = nmax - w
		f = nmax + w

		freq_centroid  = (self.freq[i:f] * self.power[i:f]).sum()
		freq_centroid /= self.power[i:f].sum()

		if (self.freq_guess == 0): self.freq_guess = freq_centroid

		return freq_centroid

	def guessParams(self):

		nmax = self.power.argmax()
		if (nmax > self.w0):
			w = self.w0
		else:
			w = nmax

		i = nmax - w
		f = nmax + w

		guess = []
		p = self.power[i:f]
		guess.append(self.freq_guess)
		freq2 = (self.freq[i:f]**2 * self.power[i:f]).sum() 
		freq2 /= self.power[i:f].sum()
		guess.append(2 * np.abs(freq2 - self.freq_guess**2)**0.5)
		guess.append(self.power[i:f].max() * 0.5 * np.pi * guess[1])
		guess.append(self.power[-100:].mean())
		guess.append(2.0)

		self.i = i
		self.f = f
		self.guess = guess

	def getRealSpectralFreq(self):

		if (self.freq_guess == 0): self.freq_guess = self.getCentroidFreq()
		if (len(self.guess) <= 4): self.guessParams()
		
		x = self.freq[self.i:self.f]
		y = self.power[self.i:self.f]

		par, cov = curve_fit(real_spectral_density, x, y, p0=self.guess[:4], maxfev=2000)
		self.guess[:par.shape[0]] = par.tolist()

		return par[0] # freq_real_spectral

	def getLorentzianFreq(self):

		if (self.freq_guess == 0): self.freq_guess = self.getCentroidFreq()
		if (len(self.guess) <= 4): self.guessParams()
		
		x = self.freq[self.i:self.f]
		y = self.power[self.i:self.f]

		par, cov = curve_fit(lorentzian, x, y, p0=self.guess[:4], maxfev=2000)
		self.guess[:par.shape[0]] = par.tolist()

		return par[0] # freq_lorentzian

	def getSoftLorentzianFreq(self):

		if (self.freq_guess == 0): self.freq_guess = self.getCentroidFreq()
		if (len(self.guess) <= 5): self.guessParams()

		x = self.freq[self.i:self.f]
		y = self.power[self.i:self.f]

		par, cov = curve_fit(soft_lorentzian, x, y, p0=self.guess[:5], maxfev=2000)
		self.guess[:par.shape[0]] = par.tolist()

		return par[0] #freq_soft_lorentz

	def getExponentialFreq(self):

		if (self.freq_guess == 0): self.freq_guess = self.getCentroidFreq()
		if (len(self.guess) <= 4): self.guessParams()

		x = self.freq[self.i:self.f]
		y = self.power[self.i:self.f]

		self.guess[2] = np.log(self.power.max())

		par, cov = curve_fit(exp_peak, x, y, p0=self.guess[:4], maxfev=2000)
		self.guess[:par.shape[0]] = par.tolist()

		return par[0] #freq_soft_lorentz

	def computePhase(self):

		# Make sure guess freq is set
		if (self.freq_guess == 0): self.freq_guess = self.getCentroidFreq()

		# Get the relative phase
		phase = np.arctan2(self.wf_im, self.wf)

		# Compute the unwrapped phase
		self.phase = np.unwrap(phase, discont=0.8*np.pi)

	def getPhaseFreq(self, poln=1):

		# Make sure guess freq is set
		if (self.phase.shape[0] == 0): self.computePhase()

		# Extract frequency with a polynomrial fit
		if (poln == 1):

			def poly(t, w, b):
				return b + w * 6.2831853 * t

		elif (poln == 2):

			def poly(t, w, b, a2):
				return b + w * 6.2831853 * t + a2 * t**2

		elif (poln >= 3):

			def poly(t, w, b, a2, a3):
				return b + w * 6.2831853 * t + a2 * t**2 + a3 * t**3

		guess = []
		guess.append(self.freq_guess)
		guess.append(self.freq_guess * (-self.fid_i))
		guess.append(0.0)
		guess.append(0.0)

		w = 10
		i = self.fid_i + w
		f = self.fid_f - w
		t = self.tm[i:f]
		p = self.phase[i:f]

		par, cov = curve_fit(poly, t, p, p0=guess[:1+poln])
		
		return par[0] # freq_phase

	def getPhaseFit(self, poln=1):

		# Make sure guess freq is set
		if (self.phase.shape[0] == 0): self.computePhase()

		# Extract frequency with a polynomrial fit
		if (poln == 1):

			def poly(t, w, b):
				return b + w * 6.2831853 * t

		elif (poln == 2):

			def poly(t, w, b, a2):
				return b + w * 6.2831853 * t + a2 * t**2

		elif (poln >= 3):

			def poly(t, w, b, a2, a3):
				return b + w * 6.2831853 * t + a2 * t**2 + a3 * t**3


		guess = []
		guess.append(self.freq_guess)
		guess.append(self.freq_guess * (-self.fid_i))
		guess.append(1.0)
		guess.append(1.0)

		w = 10
		i = self.fid_i + w
		f = self.fid_f - w
		t = self.tm[i:f]
		p = self.phase[i:f]

		par, cov = curve_fit(poly, t, p, p0=guess[:1+poln])
		
		if (poln == 1):

			fit = poly(t, par[0], par[1])

		elif (poln == 2):

			fit = poly(t, par[0], par[1], par[2])

		elif (poln >= 3):

			fit = poly(t, par[0], par[1], par[2], par[3])

		return p, fit

	def getOscillatorFreq(self):

		# Make sure guess freq is set
		if (self.freq_guess == 0): self.freq_guess = self.getCentroidFreq()

		def osc(x, f, phi, amp):
			return amp * np.sin(6.2831853 * f * x + phi)

		guess = []
		guess.append(self.freq_guess)
		guess.append(1.57) # Phase, it should handle that
		guess.append(1.0) # Amplitude should actually be close to 1.

		w = 10
		i = self.fid_i + w
		f = self.fid_f - w

		tm = self.tm[i:f]
		wv = self.wf[i:f] / self.env[i:f]

		par, cov = curve_fit(osc, tm, wv, p0 = guess)

		return par[0]

	def getFreq(self, id):

		id = self.method[id]

		if (id == 0):
			return self.getZeroCountFreq()

		elif (id == 1):
			return self.getCentroidFreq()

		elif (id == 2):
			return self.getLorentzianFreq()

		elif (id == 3):
			return self.getSoftLorentzianFreq()

		elif (id == 4):
			return self.getExponentialFreq()

		elif (id == 5):
			return self.getPhaseFreq()

		elif (id == 6):
			return self.getPhaseFreq(poln=2)

		elif (id == 7):
			return self.getPhaseFreq(poln=3)

		elif (id == 8):
			return self.getOscillatorFreq()

		else:
			print "Method ID not found.  Using 1st Order phase fit."
			return self.getPhaseFit()


	def tryAllFreqMethods(self):

		freqs = {}

		for id in self.method.keys():

			if isinstance(id, str):
				freqs[id] = self.getFreq(id)

		return freqs


"""
The goal is the test the accuracy of my FID frequency extraction methods with 
pristine waveforms generated functionally, not simulation.

"""
def ideal_fid_func(npoints, ttotal=10.0, s2n=0.01, t0=-1.0, tau2=5.0, f=23.0):

	# Set the rise time constant
	tau1 = 0.001
	omega = 2 * np.pi * f

	# Generate a random phase
	phi  = 0.0 #np.random.random_sample() * 2 * np.pi - np.pi
	# print "Using %f as phase." % phi

	# Initialize arrays
	wf = np.empty(npoints)
	tm = np.arange(t0, t0 + ttotal, float(ttotal) / float(npoints)) 

	# Create FID
	for i, t in enumerate(tm):

		if t > 0:

			wf[i] = (1 - np.exp(-t / tau1)) * np.exp(-t / tau2)
			wf[i] *= np.cos(omega * t - phi)

		if t <= 0:

			wf[i] = 0.0

	# Rescale waveform
	wf = wf / np.abs(wf).max()

	# Add noise
	wf += np.random.normal(0.0, s2n, wf.shape)

	return tm, wf



