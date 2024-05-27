#Imports
import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt

#Functions
def resample_signal(x, y, samples, high=True):
	"""
	resample_signal(x, y, samples)

	Resamples signal to be constant samples.
	
	x - the x values of the data to be resampled.
	y - the y values of the data to be resampled.
	samples - the total number of samples in the resampled signal.
	high - preserve detail but hallucinate more?

	Note: it is assumed that the resampled signal has the same x bounds as the original signal.
		The x-values of the original signal are required to be pre-sorted.

	Returns: [x_resampled, y_resampled, sampling_freq]
	"""
	x_low = x[0]
	x_high = x[-1]

	i = 0

	x_new = []
	y_new = []
	for curr_x in np.linspace(x_low,x_high,samples):
		if not high:
			if curr_x > x[i + 1]:
				i += 1
		x_new.append(curr_x)

		y_calc = y[i] + (y[i+1]-y[i])*(curr_x - x[i])/(x[i+1]-x[i])
		y_new.append(y_calc)
		if high:
			if curr_x > x[i + 1]:
				i += 1

	return [x_new, y_new, (samples-1)/(x_high - x_low)]


def butter_lowpass_filter(data, cutoff, sf, order=4): #Data to be filtered, cutoff frequency in Hz, sampling rate, polynomial order
	normalized_cutoff = 2*cutoff/sf
	b,a = butter(order, normalized_cutoff, btype='low', analog=False)
	return filtfilt(b, a, data)

def estimate_integral(x_vals, y_vals): #Estimate integral using trapezoidal rule
	"""
	area = 0
	for i in range(len(x_vals) - 1):
		area += 0.5*(y_vals[i]+y_vals[i+1])*(x_vals[i+1]-x_vals[i])
	return area
	"""
	return np.sum(0.5*(y_vals + np.roll(y_vals, 1, axis=0)) * (x_vals - np.roll(x_vals, 1, axis=0)),axis=0)

def determine_noise(molecule_wavelengths, molecule_width, xvals, yvals, xvals_new, yvals_new, samples=1000):
	l = len(xvals)

	mol_absorption = np.zeros(xvals.shape)
	for wavelength in molecule_wavelengths:
		wl = max(wavelength, 280)
		wl = min(wavelength, 1100)
		mol_absorption += np.exp(-np.power((xvals-wl)/molecule_width,2))

	mol_new = np.dstack([mol_absorption for i in range(samples)])[0]

	avgpower = estimate_integral(xvals, yvals * mol_absorption)

	power = estimate_integral(xvals_new, (yvals_new * (1 + 0.01*np.random.randn(l,samples))) * mol_new)
	diff = ((power - avgpower)/avgpower)**2
	totalNoise = np.sum(diff)
	return totalNoise*np.power(10,10)/samples


#main
if __name__ == '__main__':
	#Settings
	NUM_RESAMPLES = 10000
	LOWPASS_FREQ = 0.025
	MOL_ABSORPTION_WIDTH = 30
	INIT_WAVELENGTHS = [700,700]
	SAMPLES = 1000
	VARIANCE = 300

	#Load spectrum
	spectrum = open('sunsurf.txt', 'r')

	#Strings to float data
	xvals = []
	yvals = []
	for row in spectrum:
		cleanRow = list(map(float,row.strip().split(' ')))
		xvals.append(cleanRow[0])
		yvals.append(cleanRow[1])

	#Resample and low-pass filtering
	resampled = resample_signal(xvals, yvals, NUM_RESAMPLES, False)
	lowpass = butter_lowpass_filter(resampled[1], LOWPASS_FREQ, resampled[2])

	xvals_new = np.dstack([resampled[0] for i in range(SAMPLES)])[0]
	yvals_new = np.dstack([lowpass for i in range(SAMPLES)])[0]

	#Optimize
	print(estimate_integral(resampled[0], lowpass))

	curr = determine_noise(INIT_WAVELENGTHS,MOL_ABSORPTION_WIDTH,np.array(resampled[0]),np.array(lowpass),xvals_new,yvals_new,SAMPLES) 
	print(curr)

	final = []
	finalscore = 100000000

	l = len(INIT_WAVELENGTHS)
	for i in range(50*l):
		chngint = np.random.randint(l)
		chngvar = VARIANCE * np.random.normal()

		if i % 10 == 0:
			curr = determine_noise(INIT_WAVELENGTHS,MOL_ABSORPTION_WIDTH,np.array(resampled[0]),np.array(lowpass),xvals_new,yvals_new,SAMPLES)
		for mult in [1,-1,0.1,-1,0.1,-1,0.1,-1]:
			chngvar *= mult
			changed = True
			while changed:
				INIT_WAVELENGTHS[chngint] += chngvar
				for j in range(5):
					nc = determine_noise(INIT_WAVELENGTHS, MOL_ABSORPTION_WIDTH, np.array(resampled[0]), np.array(lowpass),xvals_new,yvals_new,SAMPLES)
					if nc >= curr-20:
						INIT_WAVELENGTHS[chngint] -= chngvar
						changed = False
						break
				if changed:
					curr = nc
		
		#curr = determine_noise(INIT_WAVELENGTHS,MOL_ABSORPTION_WIDTH,np.array(resampled[0]),np.array(lowpass),xvals_new,yvals_new,SAMPLES)

		INIT_WAVELENGTHS.sort()

		if curr < finalscore:
			finalscore = curr
			final = INIT_WAVELENGTHS

		print(INIT_WAVELENGTHS, curr, final, finalscore, i)
	print('\n')
	print(final, finalscore)

	#Plot
	plt.plot(xvals, yvals)
	plt.plot(resampled[0], lowpass)
	for wavelength in final:
		plt.axvline(x=wavelength)
	plt.show()