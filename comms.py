#Imports
import numpy as np
from numpy import cos, sin
import matplotlib.pyplot as plt

def to_rad(deg):
	return deg * np.pi/180

goldenRatio = (1 + np.sqrt(5))/2

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - 0.75*plot_radius, z_middle + 0.75*plot_radius])

#Functions
class Orbit:
	def __init__(self, a=0, e=0, i=0, lan=0, argp=0):
		self.a = a
		self.e = e
		self.i = i
		self.lan = lan
		self.argp = argp

	def genOrbit(self, minAlt, maxAlt):
		apo1 = minAlt + (maxAlt - minAlt) * np.random.random()
		apo2 = minAlt + (maxAlt - minAlt) * np.random.random()
		self.a = (apo1 + apo2)/2
		self.e = abs((apo1-apo2)/(apo1+apo2))
		self.e = min(self.e, 0.97)
		self.i = np.pi*np.random.random() - np.pi/2
		self.lan = 2*np.pi*np.random.random()
		self.argp = 2*np.pi*np.random.random()

	def calcPoint(self, v):
		d = self.a * (1-self.e * self.e)/(1+self.e*cos(v))
		x = d*cos(self.argp + v)
		y = d*sin(self.argp + v)
		y, z = y*cos(self.i), y*sin(self.i)
		x, y = x*cos(self.lan) - y*sin(self.lan), x*sin(self.lan) + y*cos(self.lan)
		return x, y, z

	def plot(self, res, color):
		X, Y, Z = self.calcPoint(np.linspace(0, 2 * np.pi, res))
		ax.plot(X, Y, Z, color)

	def mutate(self, minAlt, maxAlt):
		rands = np.random.rand(5)
		na = self.a * (0.25 * rands[0] + 0.875)
		ne = self.e * (0.25 * rands[1] + 0.875)

		ne = min(ne, (maxAlt - minAlt)/(maxAlt + minAlt))
		ne = min(ne, 0.97)
		na = max(na, minAlt/(1-ne))
		na = min(na, maxAlt/(1+ne))

		return Orbit(
			na, 
			ne,
			self.i + (0.6 * rands[2]),
			self.lan+(0.6 * rands[3]),
			self.argp+(0.6* rands[4])
		)

	def __repr__(self):
		return str([self.a, self.e, self.i, self.lan, self.argp])+"\n"

class Constellation:
	def __init__(self, orbits=[]):
		self.orbits = orbits
		self.fitness = None

	def genOrbits(self, n, minAlt, maxAlt):
		self.orbits = []
		for i in range(n):
			self.orbits.append(Orbit())
			self.orbits[-1].genOrbit(minAlt, maxAlt)

	def plot(self, res, color):
		for orbit in self.orbits:
			orbit.plot(res, color)

	def mutate(self, minAlt, maxAlt):
		return Constellation([orbit.mutate(minAlt, maxAlt) for orbit in self.orbits])

	def __repr__(self):
		out = ""
		for orbit in self.orbits:
			out += str(orbit)
		return out

	def genFibSphere(self, n):
		#http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/#:~:text=The%20Fibonacci%20lattice%20(aka%20Fibonacci,preserving%2C%20not%20distance%2Dpreserving.
		i = np.arange(0, n)
		theta = 2 * np.pi * i / goldenRatio
		phi = np.arccos(1 - 2*(i+0.5)/n)
		return cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)

	def pointsVisible(self, x1, y1, z1, x2, y2, z2):
		vx, vy, vz = x2-x1, y2-y1, z2-z1
		a = vx*vx + vy*vy + vz*vz
		b = vx*x1 + vy*y1 + vz*z1
		c = x1*x1 + y1*y1 + z1*z1 - 1

		return b*b < a*c

	def calcFitness(self, res=100, p=100):
		numS = len(self.orbits)
		dirx, diry, dirz = self.genFibSphere(p)
		self.fitness = 0
		for i in range(res):
			Ms = 2 * np.pi * np.random.rand(numS)
			E = Ms
			sX, sY, sZ = np.zeros(numS), np.zeros(numS), np.zeros(numS)
			for k in range(numS):
				for j in range(5):
					eo = self.orbits[k].e
					#print(eo)
					E[k] += (Ms[k]-(E[k]-eo*sin(E[k])))/(1-eo*cos(E[k]))
				E[k] = np.arctan2(np.sqrt(1-eo*eo)*sin(E[k]),cos(E[k])-eo)
				sX[k], sY[k], sZ[k] = self.orbits[k].calcPoint(E[k])
			graph = [[self.pointsVisible(sX[k],sY[k],sZ[k],sX[j],sY[j],sZ[j]) for k in range(numS)] for j in range(numS)]
			connected_components = [-1 for j in range(numS)]
			numcc = 0
			for j in range(numS):
				for k in range(numS):
					if graph[j][k]:
						if connected_components[k] != -1:
							connected_components[j] = connected_components[k]
				if connected_components[j] == -1:
					connected_components[j] = numcc
					numcc += 1
			numGStations = 0
			numSStations = 0
			for j in range(p):
				dx = dirx[j]
				dy = diry[j]
				dz = dirz[j]
				currS = False
				for k in range(numS):
					D = np.sqrt(sX[k]*sX[k] + sY[k]*sY[k] + sZ[k]*sZ[k])
					dot = dx * sX[k] + dy * sY[k] + dz * sZ[k]
					if dot < D * D - 1 and not currS:
						numSStations += 1
						currS = True
					if dot < -1:
						numGStations += 1
						break

			#for j in range(len(graph)):
				#for k in range(len(graph[0])):
					#if graph[j][k]:
						#ax.plot([sX[j],sX[k]],[sY[j],sY[k]],[sZ[j],sZ[k]], "black")
			self.fitness += numGStations*numSStations/(p*p*numcc)
		return self.fitness/res



c = Constellation()
c.genOrbits(4, 1.15, 10)
#c.plot(50, "red")

"""
for i in range(200):
	if i % 10 == 0:
		print('_')
	m = c.mutate(1.15, 10)
	if m.calcFitness(res=200) > c.calcFitness(res=200):
		c = m
		print(c.calcFitness(res=20))
"""

print(c)

#Number of sats
n = 4

# Make data
R = 1

u = np.linspace(0, 2 * np.pi, 20)
v = np.linspace(0, np.pi, 10)
x = R * np.outer(np.cos(u), np.sin(v))
y = R * np.outer(np.sin(u), np.sin(v))
z = R * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color=[1,0,0])
print(c.calcFitness(res=200))
#c.plot(50, "blue")

o = Orbit(10,0,0,0,0)
o.plot(50, "red")

# Set an equal aspect ratio
#ax.set_box_aspect([1,1,0.92])
set_axes_equal(ax)

plt.show()