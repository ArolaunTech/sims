"""N-Body Simulator"""
#
# This simulator predicts the motion of several particles under
# their mutual gravitation. Conversions are avoided by making all units SI.
#
"""======#END======"""

#Import libraries
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import math

#Functions
def length(x, y, z):
	return math.sqrt(x * x + y * y + z * z)

def genAcc(G, pos, mass, softening):
	#Credits to Philip Mocz.

	"""
    Calculate the acceleration on each particle due to Newton's Law 
	pos  is an N x 3 matrix of positions
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	softening is the softening length
	a is N x 3 matrix of accelerations
	"""
	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r^3 for all particle pairwise particle separations 
	inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)
	inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)

	ax = G * (dx * inv_r3) @ mass
	ay = G * (dy * inv_r3) @ mass
	az = G * (dz * inv_r3) @ mass
	
	# pack together the acceleration components
	a = np.dstack((ax,ay,az))

	return a

#Initialize Simulation
N = 100 #Number of particles
mavg = 1.0 #Average mass of a particle
mstd = 0.5 #Standard deviation of mass distribution
ivel = 10.0 #Velocity of a particle
rad = 10.0 #Distance from the origin

pos = np.random.randn(N, 3) * rad
vel = np.random.randn(N, 3)
mass = mavg + mstd * np.random.randn(N)

for i in range(N):
	l = length(vel[i][0], vel[i][1], vel[i][2]) / ivel
	vel[i][0] /= l
	vel[i][1] /= l
	vel[i][2] /= l

history = [pos]
histvel = [vel]

#Simulation constants
G = 1.0 #Newton's gravitational constant
tEnd = 10.0 #End time of simulation
dt = 0.1 #Timestep of simulation
softening = 0.05 #Softening paramater

Nt = math.ceil(tEnd / dt) #Number of timesteps

fig = plt.figure()
ax = plt.axes(xlim = (-30, 30), ylim = (-30,30))
scatter = ax.scatter(pos[:,0], pos[:,1])

for i in range(Nt): #Simulation loop
	accel = genAcc(G, pos, mass, softening)
	vel += accel[0] * dt
	pos += vel * dt
	history.append(pos)
	histvel.append(vel)


def update(i, history):
	pos = history[i]
	scatter.set_offsets(pos[:,0:2])
	return scatter,

anim = animation.FuncAnimation(fig, update, interval=10, frames = Nt, fargs=(history,))
plt.show()