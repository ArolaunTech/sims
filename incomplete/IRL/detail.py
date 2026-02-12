#"""=== High Detail Raytracing Engine using SDFs ==="""
#
#   A high detail raytracing engine.
#


"""=== Import packages ==="""
print("Importing packages (0/4)")
import numpy as np #Math and array operations
print("\x1b[1F\x1b[2KImporting packages(1/4)")
from PIL import Image #Image operations
print("\x1b[1F\x1b[2KImporting packages(2/4)")
import os.path as path
print("\x1b[1F\x1b[2KImporting packages(3/4)")
import matplotlib.pyplot as plt
print("\x1b[1F\x1b[2KImports done")

"""=== Functions & Classes ==="""
print("\nInterpreting functions")

#Utility
epsilon = 0.000001

#SDFs
def sphereSDF(center, radius):
	return lambda pos: np.linalg.norm(pos-center) - radius

def depthFromSDF(pos, direction, sdf):
	curr = sdf(pos)
	t = 0

	while curr > 0.001:
		t += curr
		curr = sdf(pos + t*direction)
		if curr > 1000:
			return -1 #No intersection
	return t #Intersection

def normalFromSDF(pos, sdf):
	out = np.array([0.0,0.0,0.0])
	for direction in [[1.0,-1.0,-1.0], [-1.0,-1.0,1.0], [-1.0,1.0,-1.0], [1.0,1.0,1.0]]:
		out += np.array(direction)*sdf(pos + epsilon*np.array(direction))
	return out/np.linalg.norm(out)

print("\x1b[1F\x1b[2KFunctions interpreted")

"""=== Define scene ==="""
print("\nDefining Scene")

#Type: list[function]
sceneSDFs = [ # SDF of various objects in scene
	sphereSDF(np.array([0,0,2]), 1),
	sphereSDF(np.array([1,0,3]), 1)
]
sceneDepths = [ # Depth for optimization
	
]
sceneBRDFsR = [] # Bi-directional reflectance function of objects in scene (red)
sceneBRDFsG = [] # Bi-directional reflectance function of objects in scene (green)
sceneBRDFsB = [] # Bi-directional reflectance function of objects in scene (blue)

#Type: list
lightdir = np.float64(np.array([1,1,-1]))
lightdir /= np.linalg.norm(lightdir)

print("\x1b[1F\x1b[2KScene defined")

"""=== After functions ==="""

def rayMarch(pos, direction): #Calculate first intersection
	depth, index = np.inf,-1
	for i, sdf in enumerate(sceneSDFs):
		curr = depthFromSDF(pos, direction, sdf)
		if curr < depth and curr != -1:
			depth = curr
			index = i
	return [depth, index]

	"""
	for i in range(1000):
		t += curr[0]
		curr = sceneSDF(pos + t*direction)
		if curr[0] < 0.001:
			return [t, curr[1]]
		if curr[0] > 1000:
			return [-1, -1]
	"""

def calcIllumination(pos, direction, maxBounces):
	depth = rayMarch(pos, direction)
	hitIndex = depth[1]
	depth = depth[0]
	col = np.dot(normalFromSDF(pos+depth*direction, sceneSDFs[hitIndex]), lightdir)
	col = np.array([col, col, col])
	return np.clip(col, 0.01, np.inf)

"""=== Image resolution ==="""
rows = 240
cols = 426

"""=== Import overwrite image ==="""
print("\nImporting Image")

inoutPath = 'detailRender.png' #Input/output file path

if path.isfile(inoutPath): #Get numpy array
	modifyImage = np.array(Image.open(inoutPath))
else:
	modifyImage = np.zeros((rows, cols, 3))
modifyImage = np.float64(modifyImage)/255

#Clear image (optional)
modifyImage = np.zeros((rows, cols, 3))

#plt.imshow(modifyImage)
#plt.show()
print("\x1b[1F\x1b[2KImage imported")

"""=== Get directions ==="""
print("\nGetting directions")

xy = np.indices((rows, cols)) #Get indices
m = min(rows, cols)

x = 2*np.float64(xy[1])/m - cols/m #Transform into coordinates
y = -2*np.float64(xy[0])/m + rows/m
z = np.ones((rows, cols))

#print(x, y)

vFOV = 60
vFOV /= 2
vFOV = np.tan(vFOV * np.pi/180)

x *= vFOV
y *= vFOV

Dd = np.sqrt(x*x + y*y + z*z) #Normalize to unit vector

x /= Dd
y /= Dd
z /= Dd

direction = np.dstack([x, y, z])

#modifyImage = direction

#plt.imshow(direction)
#plt.show()

print("\x1b[1F\x1b[2KDirections calculated")

"""=== Rendering ==="""
print('\nBeginning Rendering')

#Get row to modify
for rMod in range(rows):
	if (modifyImage[rMod][0] == 0).all():
		break
if rMod > 0:
	rMod -= 1

#Iterate over pixels and render
for r in range(rMod, rows):
	for c in range(cols):
		#For each pixel
		modifyImage[r][c] = calcIllumination(np.array([0.0,0.0,0.0]), direction[r][c], 1)
	#Save image
	img = Image.fromarray(np.uint8(np.clip(modifyImage,0,1)*255))
	img.save(inoutPath)


"""=== Save final Image ==="""
print("\x1b[1F\x1b[2KImage saved")
img = Image.fromarray(np.uint8(np.clip(modifyImage,0,1)*255))
img.save(inoutPath)