#TPT cold fusion
import matplotlib.pyplot as plt
import numpy as np
from timemodule import timer

pertimer = timer()
pertimer.tagtime() 

def getNiceness(r, c, arr):
	nice = np.zeros((r, c))
	for i in range(r):
		for j in range(c):
			if arr[i][j] == 1:
				continue
			niceness = 0
			for k in range(-1, 2):
				for l in range(-1, 2):
					ni = (i + k) % r
					nj = (j + l) % c
					if arr[ni][nj] == 0:
						continue
					niceness -= 1
					for m in range(-1, 2):
						for n in range(-1, 2):
							nni = (ni + m) % r
							nnj = (nj + n) % c
							if arr[nni][nnj] == 0:
								niceness += 1
			nice[i][j] = niceness
	return [nice, np.sum(nice)/(r*c)]

r = 4
c = 4

b = 0
bi = -1
f = 0

for i in range(2**(r*c)):
	if i in [21845, 43690, 3855, 7710, 61680, 1807, 2831, 3343, 3599, 3847, 3851, 3853, 3854, 5461, 10922, 17749, 20821, 21589, 21781, 21829, 21841, 21844, 28912, 35498, 41642, 43178, 43562, 43658, 43682, 43688, 45296, 53488, 57584, 61552, 61616, 61648, 61664, 1295, 2575, 3845, 3850, 5397, 10794, 17733, 20720, 20817, 21588, 35466, 41200, 41634, 43176, 61520, 61600]:
		continue
	patt = list(map(int, list(bin(i)[2:])))
	patt = [0 for j in range(r*c - len(patt))] + patt
	patt = np.reshape(patt, (r,c))
	g = getNiceness(r, c, patt)
	if g[1] > b:
		b = g[1]
		bi = i
		f = patt
patt = f
print(b, bi)

pattr = np.zeros((100, 100))

R = np.arange(r)
C = np.arange(c)
C,R = np.meshgrid(C, R)
print(R, C)

R2 = np.arange(100)
C2 = np.arange(100)
C2,R2 = np.meshgrid(C2, R2)
print(R2, C2)

pattr = patt[R2%r,C2%c]

plt.imshow(pattr)
print(pertimer.elapsed())
plt.show()