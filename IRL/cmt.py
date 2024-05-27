import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

delta = 6/29

def finv(t):
	if t > delta:
		return t * t * t
	return 3 * delta**2 * (t - 4/29)

def LabXYZ(L, a, b):
	#L*, a*, b* to X, Y, Z assuming standard illuminant D65.
	l = (L + 16) / 116
	return 0.950489*finv(l + a/500), finv(l), 1.0884*finv(l - b/200)

def XYZRGB(X, Y, Z):
	#X, Y, Z to R, G, B (linear)
	return 3.24*X - 1.54*Y - 0.5*Z, -0.97*X + 1.88*Y + 0.04*Z, 0.06*X - 0.2*Y + 1.06*Z

def gamma(x):
	#Linear to gamma-corrected
	if x < 0.0031308:
		return 12.92*x
	return 1.055 * x**(5/12) - 0.055

def LabRGB(Lab):
	X, Y, Z = LabXYZ(Lab[0], Lab[1], Lab[2])
	R, G, B = XYZRGB(X, Y, Z)
	R = np.clip(R, 0, 1)
	G = np.clip(G, 0, 1)
	B = np.clip(B, 0, 1)
	return [gamma(R), gamma(G), gamma(B)]

C = list(np.linspace(0, 1, num=256))

def L(x):
	return 100 * (0.2+0.75*x)

def a(x):
	return 40000 * (x*x*np.exp(-(x+1.75)**2))

def b(x):
	return 600 * (-x + 2 * x**1.9)

colores = []

for i in C:
	colores.append(LabRGB([L(i), a(i), b(i)]))

newcmp = colors.ListedColormap(colores) 

X, Y = np.meshgrid(range(500), range(100))

test = 0.002*X+0.00001*np.sin(X)*(100 - Y)**2

fig, ax = plt.subplots()

ax.imshow(list(test), vmin=0, vmax=1, cmap=newcmp)
plt.show()