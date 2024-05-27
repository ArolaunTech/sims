"""=== Custom Planetary Renderer ==="""
# Shadertoy: https://www.shadertoy.com/view/dtjfWm#
#

###Imports
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import opensimplex as osx
import sys
import time

print("Imports done")

osx.seed(7809) #Terrain seed
np.random.seed(8967) #Raytrace seed

command_args = sys.argv[1:]

###Utils
epsilon = 0.000001

###Planet Settings

##Physical parameters
ppos = np.float64(np.array([0, 0, 2]))

radius = 1.0

###Render settings
sdfFudge = 1.0/1.2 #Higher = more performance, lower = more accuracy. Set between 0 and 1, ideally close to 1.

#SDF
def planetSDF(pos):
	dcenter = np.linalg.norm(ppos-pos)
	#dnoise = 0.1*osx.noise3(ppos[0]-pos[0], ppos[1]-pos[1], ppos[2]-pos[2])
	dnoise = 0.05*np.sin((ppos[1]-pos[1])*20)
	#dnoise = min(dnoise, 0.0)

	return dcenter + dnoise - radius

def depthFromSDF(pos, direction):
	curr = planetSDF(pos)
	t = 0

	while curr > 0.001:
		t += sdfFudge*curr
		curr = planetSDF(pos + t*direction)
		if curr > 1000:
			return -1 #No intersection
	return t #Intersection

def normalFromSDF(pos):
	#From https://iquilezles.org/articles/normalsSDF/
	out = np.array([0.0,0.0,0.0])
	for direction in [[1.0,-1.0,-1.0], [-1.0,-1.0,1.0], [-1.0,1.0,-1.0], [1.0,1.0,1.0]]:
		out += np.array(direction)*planetSDF(pos + epsilon*np.array(direction))
	return out/np.linalg.norm(out)

def uniformSphere():
	p = np.array([np.random.normal(), np.random.normal(), np.random.normal()])
	return p / np.linalg.norm(p)

def calcIllumination(pos, direction, maxBounces):
	tca = np.dot(direction, ppos-pos)
	d = np.sqrt(np.dot(ppos-pos,ppos-pos)-tca*tca)
	thc = np.sqrt(atmoheight*atmoheight - d*d)

	depth = depthFromSDF(pos, direction)
	normal = normalFromSDF(pos + depth*direction)

	atmostart = 0
	atmodepth = 0

	if d < atmoheight:
		atmostart = tca-thc
		atmodepth = 2*thc
	if depth > -1:
		atmostart = min(atmostart, depth)
		atmodepth = depth - atmostart

	atmoend = atmostart + atmodepth

	atmostartp = np.float64(direction) * atmostart
	terrstartp = np.float64(direction) * depth

	terrstartalt = np.linalg.norm(ppos-terrstartp) - radius

	if depth < 0:
		terrstartalt = -1

	terrstartalt = np.array([terrstartalt, terrstartalt, terrstartalt])
	
	if depth > -1:
		if maxBounces == 0:
			return np.array([0,0,0.3])
		d = uniformSphere() + normal
		d /= np.linalg.norm(d)
		return calcIllumination(terrstartp + 0.002*normal, d, maxBounces - 1) * np.array([1.0,1.0,1.0])
		#return depthFromSDF(terrstartp + 0.002*normal, d)
	if np.dot(direction, lightdir) > 0.996:
		return np.array([100,100,100])
	return np.array([0,0,0.3])

##Atmospheric parameters
atmoheight = 0.1
atmoheight += radius

###Image settings

##Resolution
rows = 240
cols = 426

##FOV
vFOV = 60

vFOV = np.tan(vFOV*np.pi/360)

##Lighting
lightdir = np.float64(np.array([1,1,-1]))

lightdir /= np.linalg.norm(lightdir)

###Determine ray direction
xy = np.indices((rows, cols))
m = min(rows, cols)
Dx = 2*xy[1]/m-cols/m
Dy = rows/m-2*xy[0]/m
Dz = np.ones((rows, cols))

Dx *= vFOV
Dy *= vFOV

Dd = np.sqrt(Dx*Dx + Dy*Dy + Dz*Dz)

Dx /= Dd
Dy /= Dd
Dz /= Dd

direction = np.dstack((Dx, Dy, Dz))

###Render

if 'clear' in command_args:
	final = np.zeros((rows, cols, 3))
else:
	img = Image.open('planetRender.png')
	final = np.float64(np.array(img))/255

startTime = time.monotonic()
startRow = 0

for r in range(rows):
	for c in range(cols):
		if final[r][c][0] == 0:
			final[r][c] = np.array([0,0,0])

			depth = depthFromSDF(np.array([0,0,0]), direction[r][c])
			iters = 10000 if depth > -1 else 1
			for i in range(iters):
				final[r][c] += calcIllumination(np.array([0,0,0]), direction[r][c], 10)
			final[r][c] = np.clip(final[r][c]/iters, 0.01, 1)
			if 'pix' in command_args:
				img = Image.fromarray(np.uint8(255*np.clip(final, 0, 1)))
				img.save('planetRender.png')
		else:
			startRow = r
	print((time.monotonic() - startTime) * (rows - r)/(r - startRow+1))

	img = Image.fromarray(np.uint8(255*np.clip(final, 0, 1)))
	img.save('planetRender.png')

###Create image

final = np.uint8(255*np.clip(final, 0, 1))
img = Image.fromarray(final)
img.save('planetRender.png')
if 'matplot' in command_args:
	plt.imshow(final)
	plt.show()