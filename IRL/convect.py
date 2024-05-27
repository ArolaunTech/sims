"""
Convection Simulation
"""
#Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#Functions
def conserved(rho, vx, vy, P, gamma, vol):
	mass = rho * vol #Mass matrix
	mx = mass * vx #X-momentum
	my = mass * vy #Y-momentum
	E = 0.5 * mass * (vx * vx + vy * vy) + P * vol / (gamma - 1) #Energy
	return mass, mx, my, E

def primitive(mass, mx, my, E, gamma, vol):
	rho = mass / vol #Density
	vx = mx / mass #X-velocity
	vy = my / mass #Y-velocity
	P = (E - 0.5 * mass * (vx * vx + vy * vy)) * (gamma - 1) / vol
	return rho, vx, vy, P

def gradient(f, s):
	dx = (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / (2 * s)
	dy = (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2 * s)
	return dx, dy
#Initialization
s = 0.01 #Side length of a square
volume = s * s #"Volume" of a square
rows = 1
cols = 1
gamma = 5/3 #Ratio of specific heats
dt = 0.1 #Timestep

Mass = np.array([[0.0001]])
Mx = np.array([[0]])
My = np.array([[0]])
E = np.array([[0]])

t = 0 #Starting time
end = 2 #Ending time

N = (end - t) // dt

for i in range(N):
	#Main simulation loop

	#Primitive variables
	rho, vx, vy, P = primitive(Mass, Mx, My, E, gamma, volume)

	#Calculate conserved fluxes
	#Mass flux is the divergence of the momentum field
	dMass = (np.roll(Mx, 1, axis=0) - np.roll(Mx, -1, axis=0)) + (np.roll(My, 1, axis=1) - np.roll(My, -1, axis=1))
	#Mass coming in and out carries momentum and energy
	dMx = (np.roll(Mx, 1, axis=0) * np.roll(vx, 1, axis=0) - np.roll(Mx, -1, axis=0) * np.roll(vx, -1, axis=0)) + (np.roll(My, 1, axis=1) * np.roll(vx, 1, axis=1) - np.roll(My, -1, axis=1) * np.roll(vx, -1, axis=1)) + (np.roll(P, 1, axis=0) - np.roll(P, -1, axis=0))
	dMx = (np.roll(Mx, 1, axis=0) * np.roll(vy, 1, axis=0) - np.roll(Mx, -1, axis=0) * np.roll(vy, -1, axis=0)) + (np.roll(My, 1, axis=1) * np.roll(vy, 1, axis=1) - np.roll(My, -1, axis=1) * np.roll(vy, -1, axis=1)) + (np.roll(P, 1, axis=1) - np.roll(P, -1, axis=1))


	#Apply fluxes