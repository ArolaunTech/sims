"""#----- 1-D stellar simulation -----

All units are in SI to avoid unit confusion.

"""#----------------------------------

#Imports
import numpy as np
import matplotlib.pyplot as plt
from timemodule import timer
from numpy import exp

progTime = timer()
progTime.tagtime()

#Physical constants
G = 6.674e-11 #Gravitational constant
c = 299792458 #Speed of light
sigma = 5.67e-8 #Stefan-Boltzmann constant
M_sol = 1.989e+30 #Solar mass
elem = 1.60217663e-19 #Elementary charge, electronvolt

#Star settings
M = 1 #Star mass in solar masses
X = 0.71 #Hydrogen mass fraction
Z = 0.019 #Metallicity
Y = 1 - X - Z #Helium fraction

#Simulation settings
res = 10000 #Number of slices.
dt = 0.1 #Timestep
idx = np.arange(res)

"""#--- Fusion functions ---

Given a temperature in gigakelvin and abundances of its reactants, return the reaction rate in mol s^-1 cm^-3.

1 mol s^-1 cm^-3 = 1e+6 mol s^-1 m^-3

"""#------------------------

def h1_he3(T, H1): #3H-1 -> He-3
	return 4.01e-15/T**(2/3) * np.exp(-3.38/T**(1/3)) * (1+0.123*T**(1/3)+1.09*T**(2/3)+0.938*T) * H1**2

def he3he4(T, HE3): #He-3 + He-3 -> He-4 + 2H-1
	return 6.04e+10/T**(2/3) * np.exp(-12.276/T**(1/3)) * (1+0.034*T**(1/3)-0.522*T**(2/3)-0.124*T+0.353*T**(4/3)+0.213*T**(5/3)) * HE3**2

def he3be7(T, HE3, HE4): #He-3 + He-4 -> Be-7
	TA = T/(1+4.95e-2*T)
	return 5.61e+6*TA**(5/6)/T**(3/2)*np.exp(-12.826/TA**(1/3)) * HE3 * HE4

def be7spend(T, BE7, H1): #Be-7 decay into Lithium and Boron
	A = 3.11e+5/T**(2/3) * np.exp(-10.262/T**(1/3)) + 2.53e+3/T**(3/2) * np.exp(-7.306/T) * H1 * BE7 #Rate of boron creation
	B = -np.log(2)/(53.22*86400) * BE7 #Rate of beryllium decay into lithium
	return [A, B]

def dhe3(T, H1, HE3, HE4, BE7): #Rate of change in Helium-3 concentration given the concentrations of other elements.
	return he_he3(T, H1) - he3_he4(T, HE3) - he3be7(T, HE3, HE4) + be7spend(T, BE7, H1)[1]

def dbe7(T, H1, HE3, HE4, BE7): #Rate of change in Beryllium-7 concentration given the concentrations of other elements.
	spend = be7spend(T, BE7, H1)
	return he3be7(T, HE3, HE4) - spend[0] - spend[1]

def he4c12(T9, HE4): #He-4 -> C-12
	T932 = T9**(3/2)
	T913 = T9**(1/3)
	T923 = T913 * T913
	T943 = T923 * T923
	T953 = T943 * T913
	rc28 = 0.1
	if T9 <= 0.08:
		return HE4**3 * (2.90e-16*(7.40e+05/T932*exp(-1.0663/T9)+4.16e+09/T923*exp(-13.490/T913-(T9/0.098)**2)*(1.0+0.031*T913+8.009*T923+1.732*T9+49.883*T943+27.426*T953))*(1.30e+02/T932*exp(-3.3364/T9)+2.510e+07/T923*exp(-23.570/T913-(T9/0.235)**2)*(1.0+0.018*T913+5.249*T923+0.650*T9+19.176*T943+6.034*T953))*(0.01+0.2*(1.0+4.0*exp(-(0.025/T9)**3.263))/(1.0+4.0*exp(-(T9/0.025)**9.227)))+rc28*1.35e-07/T932*exp(-24.811/T9))
	return HE4** 3 * (2.79e-8/T9**3 * exp(-4.4027/T9) + rc28*1.35e-7/T932 * exp(-24.811/T9))

def c12o16(T9, C12, HE4): #C-12 -> O-16
	T913 = T9**(1/3)
	T923 = T913 * T913
	T932 = T9 * np.sqrt(T9)
	return HE4 * C12 * (1.04e+08/T9**2/(1.00+0.0489/T923)**2*exp(-32.120/T913-(T9/3.496)**2)+1.76e+08/T9**2/(1.00+0.2654/T923)**2*exp(-32.120/T913)+1.25e+03/T932*exp(-27.499/T9)+1.43e-02*T9**5*exp(-15.541/T9))

def o16ne20(T9, O16, HE4): #O-16 -> Ne-20
	T913 = T9**(1/3)
	T923 = T913 * T913
	T932 = T9 * np.sqrt(T9)
	return O16 * HE4 * (9.37e+09/T923*exp(-39.757/T913-(T9/1.586)**2)+6.21e+01/T932*exp(-10.297/T9)+5.38e+02/T932*exp(-12.226/T9)+1.30e+01*T9**2*exp(-20.093/T9))

def ne20mg24(T9, NE20, HE4): #Ne-20 -> Mg-24
	T913 = T9**(1/3)
	T923 = T913 * T913
	T943 = T923 * T923
	T953 = T943 * T913
	T912 = np.sqrt(T9)
	T932 = T9 * T912
	rc102 = 0.1
	return NE20 * HE4 * (4.11e+11/T923*exp(-46.766/T913-(T9/2.219)**2)*(1.00+0.009*T913+0.882*T923+0.055*T9+0.749*T943+0.119*T953)+5.27e+03/T932*exp(-15.869/T9)+6.51e+03*T912*exp(-16.223/T9)+rc102*4.21e+01/T932*exp(-9.115/T9)+3.20e+01/T923*exp(-9.383/T9))/(1.0+5.0*exp(-18.960/T9))

def mg24si28(T9, MG24, HE4): #Mg-24 -> Si-28
	rc121 = 0.1
	T932 = T9 * np.sqrt(T9)
	return MG24 * HE4 * (4.78e+01/T932*exp(-13.506/T9)+2.38e+03/T932*exp(-15.218/T9)+2.47e+02*T932*exp(-15.147/T9)+rc121*1.72e-09/T932*exp(-5.028/T9)+1.25e-03/T932*exp(-7.929/T9)+2.43e+01/T9*exp(-11.523/T9))/(1.0+5.0*exp(-15.882/T9))

print(mg24si28(1, 100, 50))

#Functions
def d(arr): #Differences of the array arr (must be a 1D numpy array)
	arrm = np.concatenate((arr,[0]))
	inc = np.clip(idx+1,-1,res-1)
	dec = np.clip(idx-1,-1,res-1)
	return (arrm[inc]-arrm[dec])/(inc - dec)

def calcvals(m, dm, r): #Calculates values for the next timestep based on the values of the current timestep.
	rho = 1/(4 * np.pi * r * r * d(r) / dm)
	P = np.cumsum(-G * m * rho * d(r) / (r * r))
	P -= P[-1]
	return rho, P


print(progTime.elapsed())
progTime.tagtime()

#Initialization
x = np.linspace(1/res,1,res)
m = M_sol * M * x #Mass inside of each spherical shell.
dm = np.full((res), m[0]) #Differences of m.
r = np.cbrt(x) * 696600000 * M**0.8 #Radius of shells.

rho, P = calcvals(m, dm, r)

print(P[0])

plt.plot(P)

print(progTime.elapsed())

plt.show()