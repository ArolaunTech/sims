import numpy as np
import matplotlib.pyplot as plt
import math

#Initialize simulation
rows = 100
columns = 500

simulation = np.array([ #Grid of cells. Each cell has numbers representing density, velocity (2D), temperature, pressure, and internal energy
	[[[0, [0,0], 0, 0, 0] for j in range(columns)] for i in range(rows)]
])

def sumCoeffs(r, c): #Gradient operator for position (r, c) and property p
	mult = 0
	for i in range(8):
		nr = r + dr[i]
		nc = c + dc[i]
		if nr < 0 or nr >= rows:
			continue
		if nc < 0 or nc >= columns:
			continue
		mult += coeff[i]
	return mult

def fract(x):
	return x % 1

def clamp(x, lo, hi):
	if x > hi:
		return hi
	if x < lo:
		return lo
	return x

def posclamp(r, c):
	return clamp(r, 0, rows - 1), clamp(c, 0, columns - 1)

#Constants
dr = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
dc = [-1, -1, -1, 0, 0, 0, 1, 1, 1]
coeff = [0.05, 0.2, 0.05, 0.2, 0, 0.2, 0.05, 0.2, 0.05]

#Universal gas constant
R = 8.315
#Volume of a single cell
V = 1
#Ratio of constant pressure to constant volume specific heat capacity
gamma = 1.4
#Diffusion speed
k = 1
#Viscosity
vk = 1
#Heat conductivity
hk = 1
#Molar mass
mu = 0.029
#Timestep
ts = 0.1

def stepsim():
	newsim = np.array([
		[[[0, [0, 0], 0, 0, 0] for j in range(columns)] for i in range(rows)]
	])
	for r in range(rows):
		for c in range(columns):
			#1) Diffuse density, velocity (viscous forces), and temperature (heat conduction)
			#Each cell adds its value to the other cells, which simulates diffusion
			#Heat conduction changes the energy of cells
			s = sumCoeffs(r, c)
			for i in range(8):
				nr = r + dr[i]
				nc = c + dc[i]
				if nr < 0 or nr >= rows:
					continue
				if nc < 0 or nc >= columns:
					continue
				newsim[nr][nc][0] += simulation[r][c][0] * coeff[i] / s * k
				newsim[nr][nc][1][0] += simulation[r][c][1][0] * coeff[i] / s * vk
				newsim[nr][nc][1][1] += simulation[r][c][1][1] * coeff[i] / s * vk
				newsim[nr][nc][2] += simulation[r][c][2] * coeff[i] / s * hk
			newsim[r][c][0] += simulation[r][c][0] * (1 - k)
			newsim[r][c][1][0] += simulation[r][c][1][0] * (1 - vk)
			newsim[r][c][1][1] += simulation[r][c][1][1] * (1 - vk)
			newsim[r][c][2] += simulation[r][c][2] * (1 - hk)
			#Temperature changes turn into energy changes.
			newsim[r][c][4] = simulation[r][c][4] + simulation[r][c][0] * V * R * (newsim[r][c][2] - simulation[r][c][2]) / (mu * (gamma - 1))
			#2) Calculate the pressure and temperature of the cell based on the density and internal energy.
			#Energy is transferred through heat conduction
			newsim[r][c][3] = (newsim[r][c][4] / V - 0.5 * newsim[r][c][0] * (newsim[r][c][1][0] * newsim[r][c][1][0] + newsim[r][c][1][1] * newsim[r][c][1][1]))/(1 + 1/(gamma - 1))
			newsim[r][c][2] = newsim[r][c][3] * mu / (R * newsim[r][c][0])
	#3) Simulate advection. Each cell will transfer their properties to other cells based on their velocity.
	simulation = newsim
	newsim = np.array([
		[[[0, [0, 0], 0, 0, 0] for j in range(columns)] for i in range(rows)]
	])
	for r in range(rows):
		for c in range(columns):
			#For each cell, transfer attributes.
			nr = r + simulation[r][c][1][0] * ts
			nc = c + simulation[r][c][1][1] * ts
			nr, nc = posclamp(nr, nc)

	simulation = newsim