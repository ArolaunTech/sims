import numpy as np
import noise
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from PIL import Image
from timemodule import timer
from scipy import ndimage as ndim

pertimer = timer()
pertimer.tagtime()

def bias(x, a):
	if a == 0:
		return x
	return (np.exp(a*x)-1)/(np.exp(a)-1)

res = 1000
low = [1,0]
upp = [3,2]
X = np.linspace(low[0],upp[0],res)
Y = np.linspace(low[1],upp[1],res)
Rx = np.linspace(-5000,5000,res)
Ry = np.linspace(-5000,5000,res)
A = np.arange(0,1000,1)
A, A2 = np.int64(np.meshgrid(A, A))
La = np.linspace(-90,90,res)
Lo, La = np.meshgrid(La, -La)
Rx, Ry = np.meshgrid(Rx, Ry)
X, Y = np.meshgrid(X,Y)

Z = np.zeros(X.shape)
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
		#self.u = np.where(D == 0, 0, strength/(2*np.pi) * (self.C-c)/D)
		self.u = np.zeros(self.C.shape)
		self.v = np.where(D == 0, 0, -strength/(2*np.pi) * (self.R-r)/D)
		self.psi = np.where(D == 0, 0, strength/(4*np.pi) * np.log(D))
		return self

nf = np.vectorize(noise.snoise2)

mountain = (0.5*nf(0.5*X+9,0.5*Y+9)+0.5)*(1-abs(nf(0.5*X,0.5*Y)))

terrain = mountain
multh = 0.5
mults = 0.7
for i in range(9):
	terrain += multh * nf(mults*X, mults*Y)
	multh *= 0.5
	mults *= 2

terrain -= 0.4

terrain = bias(terrain/terrain.max(),2.5)*terrain.max()

a = max(terrain.max(),-(terrain.min()))

sediment = 0.01 * np.ones(X.shape)
rock = terrain - sediment
water = 0 * sediment

terrain = np.where(terrain < 0, 0, terrain)

dFdx = (np.roll(terrain, -1, axis=1)-np.roll(terrain, 1, axis=1))*256/2
dFdy = (np.roll(terrain, -1, axis=0)-np.roll(terrain, 1, axis=0))*256/2

R = -dFdx
G = dFdy
B = np.ones(X.shape)

mag = np.sqrt(R*R+G*G+B*B)

print(str(pertimer.elapsed())+'s preparation')
pertimer.tagtime()

partslope = terrain.copy()
terrain = rock + sediment
recfrac = 0.03
numiters = 500000
recnum = int(numiters * (1 - recfrac))
recwater = np.zeros(X.shape)

for i in range(numiters):
	#New rain
	rainX = np.random.randint(res)
	rainY = np.random.randint(res)
	curr = 0

	for j in range(1000):
		water[rainX][rainY] += 1

		recwater[rainX][rainY] += 0.9996**(numiters-i)
		partslope[rainX][rainY] += 1/60000

		w = partslope[rainX-1][rainY]
		e = partslope[(rainX+1)%res][rainY]
		s = partslope[rainX][rainY-1]
		n = partslope[rainX][(rainY+1)%res]
		c = partslope[rainX][rainY]

		erode = min(0.001, sediment[rainX][rainY], max(0, c - min(w,e,s,n))) - 0.5 * curr/mag[rainX][rainY]
		terrain[rainX][rainY] -= erode
		sediment[rainX][rainY] -= erode
		partslope[rainX][rainY] -= erode
		curr += erode

		c = partslope[rainX][rainY]
		if w < min(e, s, n, c):
			rainX -= 1
			continue
		if e < min(w, s, n, c):
			rainX += 1
			rainX %= res
			continue
		if s < min(w, e, n, c):
			rainY -= 1
			continue
		if n < min(w, e, s, c):
			rainY += 1
			rainY %= res
			continue
		break
terrain -= np.minimum(sediment, water / 60000)
sediment -= np.minimum(sediment, water / 60000)

def slat(L, wsi): #Season transformation
	return np.sign(L)*np.maximum(np.zeros(X.shape), np.abs(L)-wsi)

temp = np.cos(La*np.pi/180)
recwaterthreshold = 0.5
waterarea = np.logical_or(terrain < 0, recwater > recwaterthreshold)
temp = np.where(waterarea, 0.4 + 0.6 * temp, 0.2 + 0.8 * temp)
temp = -45 + 65 * temp
temp = np.where(waterarea, temp - 3, temp - 90 * terrain)
ntempfunc = 15*(1 + 0.1*nf(97+X,97+Y) + 0.05*nf(109+5*X,109+5*Y) + 0.025*nf(39+25*X,39+25*Y))
temp += ntempfunc

wintemp = np.cos(slat(La, -10)*np.pi/180)
wintemp = np.where(waterarea, 0.4 + 0.6 * wintemp, 0.2 + 0.8 * wintemp)
wintemp = -45 + 65 * wintemp
wintemp = np.where(waterarea, wintemp - 3, wintemp - 90 * terrain)
wintemp += ntempfunc

sumtemp = np.cos(slat(La, 10)*np.pi/180)
sumtemp = np.where(waterarea, 0.4 + 0.6 * sumtemp, 0.2 + 0.8 * sumtemp)
sumtemp = -45 + 65 * sumtemp
sumtemp = np.where(waterarea, sumtemp - 3, sumtemp - 90 * terrain)
sumtemp += ntempfunc

print(np.max(temp), np.min(temp))

print(str(pertimer.elapsed())+'s terrain and temp')
pertimer.tagtime()

l = LagrangeSolution(A, A2)
t = np.maximum(terrain,0)/np.max(terrain)
for i in range(250):
	r = np.random.randint(res)
	c = np.random.randint(res)
	if terrain[r][c] <= 0:
		continue

	for j in range(100):
		w = terrain[r-1][c]
		e = terrain[(r+1)%res][c]
		s = terrain[r][c-1]
		n = terrain[r][(c+1)%res]
		c2 = terrain[r][c]

		if w >= max(e, s, n, c2):
			r -= 1
			continue
		if e >= max(w, s, n, c2):
			r += 1
			r %= res
			continue
		if s >= max(w, e, n, c2):
			c -= 1
			continue
		if n >= max(w, e, s, c2):
			c += 1
			c %= res
			continue
		break

	l += LagrangeSolution(A, A2).source(t[r][c]*10,r,c)

#Find wind
rotratemult = 0.5
wsm = -0.05
a = 0.5
b = 0.5
ws = b + (1-b) * np.cos(np.abs(np.pi*La/180)-a)

def h(x):
	return 60 * np.floor(x/60) + 30

u = rotratemult*ws*(np.cos(h(La)*np.pi/180)-np.cos(np.pi*La/180))

def g(x):
	return np.where(x < 90, 1 - ((x-30) % 60)/30, -1)

v = ws * wsm * g(La)

u += l.u
v += l.v

print(str(pertimer.elapsed())+'s wind')
pertimer.tagtime()

#Water
velocity = np.sqrt(u*u + v*v)

windmult = 6

newvel = np.minimum(1/3,velocity*windmult)

un = u * newvel/velocity
vn = v * newvel/velocity

uvs = [-1,0,1,-1,0,1,-1,0,1]
vvs = [-1,-1,-1,0,0,0,1,1,1]
weights = [1/36,1/9,1/36,1/9,4/9,1/9,1/36,1/9,1/36]
dots = [un*uvs[i] + vn*vvs[i] for i in range(9)]
eq = [weights[i]*(1+3*dots[i]+4.5*dots[i]*dots[i]+1.5*velocity**4) for i in range(9)]
seq = sum(eq)
for i in range(9):
	eq[i]/=seq

theta = 25 + 190 * velocity
tkelvin = temp + 273.15
pws = np.exp(77.345+0.0057*tkelvin-7235/tkelvin)/np.power(tkelvin,8.2)
pa = 101325
maxh = 0.62198*pws/(pa-pws)
humidity = np.zeros(X.shape) #Humidity in atmosphere
dt = 0.01
chang = np.minimum(pws/np.max(maxh),theta*dt)
chang = np.where(waterarea,chang,chang*0.3)
rainmultpress = (1+0.5*np.cos(La*np.pi/30))
rainmulthum = 1 - 10/pws

def ceilchang(x):
	return np.sign(x) * np.ceil(np.abs(x))

def lerp(x, a, b):
	return b*x + a*(1-x)

print(str(pertimer.elapsed())+'s water - humidity loop')
pertimer.tagtime()

for i in range(250):
	humidity += chang*(maxh-(humidity/(pa-humidity)))
	humidity *= rainmulthum
	new = np.zeros(X.shape)
	for j in range(9):
		new += eq[j] * np.roll(humidity,(uvs[j],vvs[j]),axis=(1,0))
	humidity = new

humidity /= pws
humidity *= rainmultpress
humidity /= np.max(np.where(waterarea,0,humidity))
#humidity *= 10

waterblurred = ndim.gaussian_filter(np.where(waterarea, np.ones(X.shape), np.zeros(X.shape)), sigma = 10)

print(str(pertimer.elapsed())+'s humidity loop')
pertimer.tagtime()

#Biomes
biomes = np.int64(29*np.ones(X.shape))
biomes = np.where(waterarea, 28*np.ones(X.shape), biomes)

#Rainforest
biomes = np.where(
	np.logical_and(wintemp>18, np.logical_and(biomes == 29, humidity>0.33)),
	np.zeros(X.shape),
	biomes
)

biomes = np.where(
	np.logical_and(np.logical_and(wintemp>18, sumtemp < 25), np.logical_and(biomes == 29, humidity>0.25)),
	np.zeros(X.shape),
	biomes
)

#Monsoon
biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, wintemp > 25), humidity>0.25),
	np.ones(X.shape),
	biomes
)

#Desert
biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, wintemp > 10), humidity<0.061),
	3*np.ones(X.shape),
	biomes
)

#Savanna
biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, wintemp > 25), humidity>0.1),
	2*np.ones(X.shape),
	biomes
)

#Steppe
biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, wintemp > 5), np.logical_and(humidity>0.061, humidity<0.1)),
	5*np.ones(X.shape),
	biomes
)

#Ice cap
biomes = np.where(
	np.logical_and(biomes==29, wintemp < -20),
	27*np.ones(X.shape),
	biomes
)

#Tundra
biomes = np.where(
	np.logical_and(biomes==29, wintemp < -10),
	26*np.ones(X.shape),
	biomes
)

#Taiga
biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, wintemp < 5), humidity > 0.05),
	24*np.ones(X.shape),
	biomes
)

#Warm tundra
biomes = np.where(
	np.logical_and(biomes==29, sumtemp < 5),
	26*np.ones(X.shape),
	biomes
)

#Rainforest
biomes = np.where(
	np.logical_and(biomes == 29, humidity > 0.25),
	25*np.ones(X.shape),
	biomes
)

#Mediterranean
biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, np.logical_and(wintemp > 22,wintemp < 30)), np.logical_and(humidity > 0.1,np.logical_and(humidity<0.15, waterblurred > 0.001))),
	7*np.ones(X.shape),
	biomes
)

biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, np.logical_and(wintemp > 15,wintemp < 20)), np.logical_and(humidity> 0.2, waterblurred > 0.001)),
	8*np.ones(X.shape),
	biomes
)

#Grassland
biomes = np.where(
	np.logical_and(biomes == 29, humidity < 0.15),
	10*np.ones(X.shape),
	biomes
)
"""
#Oceanic
biomes = np.where(
	np.logical_and(np.logical_and(biomes==29, np.logical_and(wintemp > 5,wintemp < 20)), np.logical_and(humidity> 0.2, waterblurred > 0.0001)),
	12*np.ones(X.shape),
	biomes
)
"""
#Forest
biomes = np.where(biomes == 29, 11 * np.ones(X.shape), biomes)

print(str(pertimer.elapsed())+'s biomes')
pertimer.tagtime()

terrain = np.where(terrain < 0, 0, terrain)

dFdx = (np.roll(terrain, -1, axis=1)-np.roll(terrain, 1, axis=1))*256/2
dFdy = (np.roll(terrain, -1, axis=0)-np.roll(terrain, 1, axis=0))*256/2

occpx = 1/(1+10000*(np.roll(terrain,-1,axis=1)-terrain)**2)
occmx = 1/(1+10000*(np.roll(terrain,1,axis=1)-terrain)**2)
occpx = 1/(1+10000*(np.roll(terrain,-1,axis=0)-terrain)**2)
occmx = 1/(1+10000*(np.roll(terrain,1,axis=0)-terrain)**2)

occ = 0.25 * (occpx + occmx + occpx + occmx)

R = -dFdx
G = dFdy
B = np.ones(X.shape)

mag = np.sqrt(R*R+G*G+B*B)

R /= mag
G /= mag
B /= mag

ldir = [-0.577,0.577,0.577]
light = R*ldir[0]+G*ldir[1]+B*ldir[2]

ambient = 0.2 * occ

light = np.where(light<ambient,ambient,light)

curr = cm.get_cmap("terrain").copy()
curr.set_bad(color="blue")

light = np.where(recwater > recwaterthreshold, ldir[2], light)

maxval = np.max(terrain)
minval = -0.3 * maxval

one = np.ones(X.shape)
zero = np.zeros(X.shape)

terrain = rock + sediment
tsave = (terrain - np.min(terrain))/(np.max(terrain) - np.min(terrain))
print((-np.min(terrain))/(np.max(terrain)-np.min(terrain)))
terrain = np.where(terrain < 0, 0, terrain)

R = {
	0: np.array([20, 82, 45, 255]),
	1: np.array([86, 124, 40, 255]),
	2: np.array([152, 166, 34, 255]),
	3: np.array([171, 143, 107, 255]),
	5: np.array([130, 130, 70, 255]),
	7: np.array([50, 50, 25, 255]),
	8: np.array([50, 50, 25, 255]),
	10: np.array([70, 166, 34, 255]),
	11: np.array([30, 123, 68, 255]),
	24: np.array([10, 42, 23, 255]),
	25: np.array([20, 60, 30, 255]),
	26: np.array([105, 105, 90, 255]),
	27: np.array([240, 240, 240, 255]),
	28: np.array([0, 0, 153, 255])
}

coloresdebiomes = np.zeros((1000, 1000, 4))
for r in range(1000):
	for c in range(1000):
		coloresdebiomes[r][c] = R[biomes[r][c]]

color = coloresdebiomes / 255
snow = np.logical_and(temp < 0, mag < 2)
color = np.where(np.dstack((waterarea, waterarea, waterarea, waterarea)),np.dstack((zero, zero, one, one)),color)
#color = np.where(np.dstack((recwater, recwater, recwater, recwater)) > recwaterthreshold, np.dstack((zero, zero, one, one)), color)
color = np.where(np.dstack((snow,snow,snow,snow)), np.dstack((1.7*one, 1.7*one, 1.7*one, one)),color)
im = Image.fromarray(np.uint8(255*np.minimum(np.dstack((light,light,light,one))*color,np.dstack((one,one,one,one)))))
im.save("shade.png")

filtamounts = 5
ndim.gaussian_filter(coloresdebiomes[:,:,0], filtamounts, output=coloresdebiomes[:,:,0])
ndim.gaussian_filter(coloresdebiomes[:,:,1], filtamounts, output=coloresdebiomes[:,:,1])
ndim.gaussian_filter(coloresdebiomes[:,:,2], filtamounts, output=coloresdebiomes[:,:,2])
ndim.gaussian_filter(coloresdebiomes[:,:,3], filtamounts, output=coloresdebiomes[:,:,3])

color = coloresdebiomes / 255
snow = np.logical_and(temp < 0, mag < 2)
color = np.where(np.dstack((waterarea, waterarea, waterarea, waterarea)),np.dstack((zero, zero, one, one)),color)
#color = np.where(np.dstack((recwater, recwater, recwater, recwater)) > recwaterthreshold, np.dstack((zero, zero, one, one)), color)
color = np.where(np.dstack((snow,snow,snow,snow)), np.dstack((1.7*one, 1.7*one, 1.7*one, one)),color)
im = Image.fromarray(np.uint8(255*np.minimum(np.dstack((light,light,light,one))*color,np.dstack((one,one,one,one)))))
im.save("shade2.png")

im2 = Image.fromarray(np.uint8(255 * np.dstack((tsave, tsave, tsave, one))))
im2.save("elevation.png")

tsave = np.where(recwater > recwaterthreshold, 0, tsave)
im3 = Image.fromarray(np.uint8(255 * np.dstack((tsave, tsave, tsave, one))))
im3.save("elevlake.png")

temp -= temp.min()
temp /= temp.max()
im3 = Image.fromarray(np.uint8(255 * np.dstack((temp, temp, temp, one))))
im3.save("temperature.png")

im3 = Image.fromarray(np.uint8(255 * np.dstack((0*humidity, 0.5*humidity/humidity.max(), humidity/humidity.max(), one))))
im3.save("humidity.png")

waterplot = np.where(waterarea,np.nan,0)
terrain = np.where(terrain==0, np.nan, terrain)

plot = biomes

spacing = 25
plt.quiver(Rx[0:res:spacing,0:res:spacing], Ry[0:res:spacing,0:res:spacing], u[0:res:spacing,0:res:spacing], v[0:res:spacing,0:res:spacing],scale = 1)
plt.imshow(np.int64(coloresdebiomes), cmap=curr, extent=[-5000,5000,-5000,5000])
plt.colorbar()
print(str(pertimer.elapsed())+'s plot')
plt.show()