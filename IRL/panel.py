import numpy as np
import matplotlib.pyplot as plt
import math
import timemodule
from scipy import integrate, linalg

times = timemodule.timer()
times.tagtime()

rows = 200
cols = 200
R, C = np.meshgrid(np.linspace(-5,5,cols), np.linspace(-5,5,rows))
Z = 0*R

class LagrangeSolution:
	def __init__(self, R, C, u=Z, v=Z, psi=Z):
		self.u = u
		self.v = v
		self.psi = psi
		self.R = R
		self.C = C

	def __add__(self, other):
		return LagrangeSolution(self.R, self.C, self.u+other.u, self.v+other.v, self.psi+other.psi)

	def source(self, strength, r, c):
		"""
		LagrangeSolution.source(strength, r, c)
	
		Turns the LagrangeSolution into a source/sink.

		------------
		Inputs
		------------

		strength - the strength of the source. Input negative values for a sink.

		r, c - the position of the source.

		------------
		Returns
		------------

		None
		"""
		s_mod = strength/(2*np.pi)
		D = (self.R-r)**2 + (self.C-c)**2
		self.u = np.where(D==0,0,s_mod*(self.R-r)/D)
		self.v = np.where(D==0,0,s_mod*(self.C-c)/D)
		self.psi = s_mod * np.arctan2(self.C-c, self.R-r)
		return self

	def freestream(self, r_strength, c_strength):
		"""
		LagrangeSolution.freestream(r_strength, c_strength)

		Turns the LagrangeSolution into a freestream.

		------------
		Inputs
		------------

		r_strength - the x-velocity of the freestream.

		c_strength - the y-velocity of the freestream.

		------------
		Returns
		------------

		None
		"""
		zero = 0*self.R
		self.u = zero + r_strength
		self.v = zero + c_strength
		self.psi = self.C * r_strength - self.R * c_strength
		return self

	def doublet(self, strength, r, c, theta=0):
		"""
		LagrangeSolution.doublet(strength, r, c, theta)

		Turns the LagrangeSolution into a doublet.

		------------
		Inputs
		------------

		strength - The strength of the doublet.

		r, c - The position of the doublet.
		
		theta (optional, defaults to 0) - The rotation of the doublet.

		------------
		Returns
		------------

		None
		"""
		dr = self.R-r
		dc = self.C-c
		drp = dr*np.cos(theta/2) + dc*np.sin(theta/2)
		dcp = dc*np.cos(theta/2) - dr*np.sin(theta/2)
		s = -strength/(2*np.pi)
		D = drp**2 +dcp**2
		self.u = np.where(D==0, 0, s*(drp**2-dcp**2)/D**2)
		self.v = np.where(D==0, 0, 2*s*drp*dcp/(D*D))
		self.psi = np.where(D==0, 0, s*dcp/D)
		return self

	def vortex(self, strength, r, c):
		'''
		LagrangeSolution.vortex(strength, r, c)

		Turns the LagrangeSolution into a vortex.

		------------
		Inputs
		------------

		strength - The strength of the vortex.

		r, c - The position of the vortex.

		------------
		Returns
		------------

		None
		'''
		D = (self.R-r)**2 + (self.C-c)**2
		self.u = strength/(2*np.pi) * (self.C-c)/D
		self.v = -strength/(2*np.pi) * (self.R-r)/D
		self.psi = strength/(4*np.pi) * np.log(D)
		return self

vinf = 0.2
alpha = 0
l = LagrangeSolution(R, C).doublet(3,0,0) + LagrangeSolution(R, C).freestream(vinf*np.cos(alpha), vinf*np.sin(alpha))
u, v, psi = l.u, l.v, l.psi

fig, ax = plt.subplots(figsize=(5,5))
plt.streamplot(R, C, u, v, density=2, linewidth=1, arrowsize=2, arrowstyle='->')
plt.contour(R, C, psi, levels=[0], colors='#CD2305')
plt.title("Fluid flow")
plt.grid()
plt.imshow(1-(u**2+v**2)/vinf**2, vmin=-1, vmax=1, cmap='inferno', extent=[-5,5,-5,5])
#plt.imshow(psi, vmin=-1, vmax=1)

print(times.elapsed())

plt.show()