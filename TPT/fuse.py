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

r = 5
c = 4

b = 0
bi = -1
f = 0

for i in range(2**(r*c)):
        patt = list(map(int, list(bin(i)[2:])))
        patt = [0 for j in range(r*c - len(patt))] + patt
        patt = np.reshape(patt, (r,c))
        g = getNiceness(r, c, patt)
        if g[1] > b:
                b = g[1]
                bi = i
                f = patt
        if i % 1000 == 0:
                print(i)
	
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
