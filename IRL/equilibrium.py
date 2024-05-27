import numpy
import math

rho = 1
E = 1
Momx = 1
Momy = 0

out = [
	0, 0, 0,
	0, 0, 0,
	0, 0, 0
]

def clamp(x, lo, hi):
	if x < lo:
		return lo
	if x > hi:
		return hi
	return x

E /= rho

E = clamp(E, 0, 2)

po = E - 0.5 * E * E
pd = 0.25 * E * E
pc = 1 - E + 0.25 * E * E

po *= rho
pd *= rho
pc *= rho

ca = math.atan(pd / (pd + po))

def maxRad(t):
	mt = (t + ca) % (0.5 * math.pi) - ca
	if mt < ca:
		return (pd + po) / math.cos(-ca + (t + ca) % (0.5 * math.pi))
	else:
		return math.sqrt(2) * (pd + 0.5 * po) / math.cos(-0.25 * math.pi + (t) % (0.5 * math.pi))

maxr = maxRad(math.atan2(Momy, Momx))
if math.sqrt(Momx * Momx + Momy * Momy) > maxr:
	maxr /= math.sqrt(Momx * Momx + Momy * Momy)
	Momx *= maxr
	Momy *= maxr

out[4] = pc

def c(l, r): #Returns closest element to 0 within range [l, r]
	if r < 0:
		return r
	if l > 0:
		return l
	return 0

Mxo = c(Momx - pd, Momx + pd)
Myo = c(Momy - pd, Momy + pd)

ro = 0.25 * (po - abs(Mxo) - abs(Myo))

out[1] += ro
out[3] += ro
out[5] += ro
out[7] += ro

if Mxo < 0:
	out[3] -= Mxo
else:
	out[5] += Mxo

if Myo < 0:
	out[7] -= Myo
else:
	out[1] += Myo

Mxo = Momx - Mxo
Myo = Momy - Myo

dx = 0.5 * (Mxo - Myo)
dy = 0.5 * (Mxo + Myo)

ro = 0.25 * (pd - abs(dx) - abs(dy))

out[0] += ro
out[2] += ro
out[8] += ro
out[6] += ro

if dx < 0:
	out[0] -= dx
else:
	out[8] += dx
if dy < 0:
	out[6] -= dy
else:
	out[2] += dy

print(out)