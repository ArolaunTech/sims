"""
Implementation of the Lattice Boltzmann method within python.

All units are SI to avoid unit conversions.
"""
#Import libraries

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

#Fluid constants
mol_m = 0.029 #Molar mass in kilograms per mole

#Physical constants
R = 8.315 #Universal gas constant in Joules per mole per Kelvin

#Simulation constants
s = 0.5 #Side length of cell in meters

V = s * s #Volume of cell in square meters

rows = 5 #Size of the simulation area
cols = 5

#Functions

def newsim():
	return [[[0, 0, [0, 0, 0, 0, 0, 0, 0, 0, 0]] for j in range(cols)] for i in range(rows)]

def clamp(x, lo, hi):
	if x > hi:
		return hi
	if x < lo:
		return lo
	return x

#Initialize simulation

rho = 1.225 #Average density in kilograms per square meter

vel = [ #Average velocity
	0, 0, 0,
	0, 1, 0,
	0, 0, 0
]
movy = [
	-1, 0, 1,
	-1, 0, 1,
	-1, 0, 1
]
movx = [
	-1, -1, -1,
	0, 0, 0,
	1, 1, 1
]
kmatrix = [
	2, 1, 2,
	1, 0, 1,
	2, 1, 2
]
svel = sum(vel)
for i in range(9):
	vel[i] /= svel

E = 1 #Average energy

sim = [[[rho * V, E, vel] for j in range(cols)] for i in range(rows)]

#Simulate
nsteps = 10 #Number of simulation steps
ts = 0.1 #Timestep in seconds
df = 0.5
visc = 0.5
cf = 0.01
ci = 100

diffuse = [
	[5,7,4,8],
	[3,7,4,6],
	[1,3,4,0],
	[1,5,4,2],
	[0,4,1,3],
	[2,4,1,5],
	[6,4,3,7],
	[8,4,5,7]
]

history = []
print(sim[0][0])

for n in range(nsteps):
	"""Streaming"""
	new = newsim()
	for r in range(rows):
		for c in range(cols):
			for vindex in range(9):
				nr = r + movx[vindex]
				nc = c + movy[vindex]
				nr = clamp(nr, 0, rows - 1)
				nc = clamp(nc, 0, cols - 1)
				aindex = vindex
				if nr == r and nc == c:
					aindex = 8 - vindex
				new[nr][nc][2][aindex] += sim[r][c][0] * sim[r][c][2][vindex]
				new[nr][nc][0] += sim[r][c][0] * sim[r][c][2][vindex]
				new[nr][nc][1] += sim[r][c][1] * sim[r][c][2][vindex]
	for r in range(rows):
		for c in range(cols):
			svel = sum(new[r][c][2])
			for vindex in range(9):
				new[r][c][2][vindex] /= svel
	"""Viscosity and Diffusion of density"""
	sim = new
	new = newsim()
	for r in range(rows):
		for c in range(cols):
			new[r][c][0] += (1 - df) * sim[r][c][0]
			new[r][c][1] += (1 - df) * sim[r][c][1]
			for vindex in range(9):
				new[r][c][2][vindex] += (1 - visc) * sim[r][c][2][vindex]
				if vindex == 4:
					continue
				nr = r + movx[vindex]
				nc = c + movy[vindex]
				nr = clamp(nr, 0, rows - 1)
				nc = clamp(nc, 0, cols - 1)
				new[nr][nc][0] += df * sim[r][c][0] / 8
				new[nr][nc][1] += df * sim[r][c][1] / 8
				for vindex2 in range(9):
					new[nr][nc][2][vindex2] += visc * sim[r][c][2][vindex2] / 8
	for r in range(rows):
		for c in range(cols):
			svel = sum(new[r][c][2])
			for vindex in range(9):
				new[r][c][2][vindex] /= svel
	sim = new
	"""Collision"""
	for r in range(rows):
		for c in range(cols):
			KE = 0
			for vindex in range(9):
				KE += kmatrix[vindex] * sim[r][c][2][vindex] * sim[r][c][0]
			T = sim[r][c][1] - KE
			newvels = sim[r][c][2].copy()
			for vindex in range(9):
				newvels[vindex] += T

			"""Debug"""
			if r == 0 and c == 0:
				print(T)
				print(sim[0][0][2])

			for vindex in range(ci):
				newnewvels = newvels.copy()
				for reac in diffuse:
					am = newvels[reac[0]] * newvels[reac[1]]
					newnewvels[reac[0]] -= am * cf
					newnewvels[reac[1]] -= am * cf
					newnewvels[reac[2]] += am * cf
					newnewvels[reac[3]] += am * cf
				newvels = newnewvels.copy()

			if r == 0 and c == 0:
				print(newvels)

			svel = sum(newvels)
			for vindex in range(4):
				mi = min(newvels[i], newvels[8 - i], T)
				newvels[i] -= mi
				newvels[8 - i] -= mi
			svel = sum(newvels)
			newvels[4] = max(0, newvels[4] - svel + 1)

			svel = sum(newvels)
			for vindex in range(9):
				newvels[vindex] /= svel

			sim[r][c][2] = newvels.copy()
	history.append(sim)


def f(i, r, c):
	KE = 0
	for vindex in range(9):
		KE += kmatrix[vindex] * sim[r][c][2][vindex] * sim[r][c][0]
	return history[i][r][c][1] - KE

fig,ax = plt.subplots()
ims = []
for i in range(len(history)):
	im = ax.imshow([[f(i, r, c) for c in range(cols)] for r in range(rows)], animated=True, vmin=0, vmax=1)
	if i == 0:
		ax.imshow([[f(i, r, c) for c in range(cols)] for r in range(rows)], vmin=0, vmax=1)
	ims.append([im])
ani = anim.ArtistAnimation(fig, ims, interval = 500, blit = True, repeat_delay = 1000)
plt.show()