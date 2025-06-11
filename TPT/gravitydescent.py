from copy import deepcopy
import random
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.linalg import lstsq
from datetime import datetime

ROWS = 96
COLS = 153

def calcDiff(a, b):
	out = 0
	for i in range(len(a)):
		for j in range(len(a[0])):
			out += (a[i][j][0] - b[i][j][0]) ** 2 + (a[i][j][1] - b[i][j][1]) ** 2
	return out

target = [[[0,0] for j in range(COLS)] for i in range(ROWS)]
include = [[1 for j in range(COLS)] for i in range(ROWS)]
for r in range(ROWS):
	for c in range(COLS):
		if (r-48)**2 + (c-76)**2 < 225:
			target[r][c][0] = -1 + 0.02*(48-r)
			target[r][c][1] = 0 + 0.02*(76-c)

			include[r][c] = 0

kernel = [[[0,0] for j in range(2 * COLS + 1)] for i in range(2 * ROWS + 1)]
for r in range(2 * ROWS + 1):
	for c in range(2 * COLS + 1):
		x = r - ROWS - 1
		y = c - COLS - 1
		if (x == 0) and (y == 0):
			continue
		dist = x**2 + y**2
		dist = dist * math.sqrt(dist)
		kernel[r][c][0] = -34.17 * x/dist
		kernel[r][c][1] = -34.17 * y/dist

targetMatrix = np.array(
	[[target[r][c][0]] for c in range(COLS) for r in range(ROWS)] +
	[[target[r][c][1]] for c in range(COLS) for r in range(ROWS)]
)

gravMatrix = np.array(
	[[(kernel[r-sr+ROWS+1][c-sc+COLS+1][0] if include[sr][sc] else 0) for sc in range(COLS) for sr in range(ROWS)] for c in range(COLS) for r in range(ROWS)] +
	[[(kernel[r-sr+ROWS+1][c-sc+COLS+1][1] if include[sr][sc] else 0) for sc in range(COLS) for sr in range(ROWS)] for c in range(COLS) for r in range(ROWS)]
)

solutionMatrix, res, rnk, s = lstsq(gravMatrix, targetMatrix, lapack_driver='gelsy', check_finite=False)
solution = [[solutionMatrix[c + r * COLS][0]*256 for c in range(COLS)] for r in range(ROWS)]

outputFile = ("gravityMap " + str(datetime.now()).replace(".", "_") + ".txt").replace(" ", "_")
f = open(outputFile, "w")
f.write(str(ROWS))
f.write(" ")
f.write(str(COLS))
f.write("\n")
for r in range(ROWS):
	for c in range(COLS):
		f.write(str(solution[r][c]))
		f.write("\n")

f.close()

minimum = min(map(min, solution))
maximum = max(map(max, solution))

dist = max(abs(minimum), abs(maximum))

fig, ax = plt.subplots()

plt.imshow(solution, cmap="RdBu", vmin=-dist, vmax=dist)
plt.colorbar()
plt.title("Mass distribution for given gravity distribution")
plt.xlabel("Columns")
plt.ylabel("Rows")
spacing = 2
ax.quiver(
	[[c for c in range(0, COLS, spacing)] for r in range(0, ROWS, spacing)],
	[[r for c in range(0, COLS, spacing)] for r in range(0, ROWS, spacing)],
	[[target[r][c][1]+1e-6 for c in range(0, COLS, spacing)] for r in range(0, ROWS, spacing)],
	[[target[r][c][0]+1e-6 for c in range(0, COLS, spacing)] for r in range(0, ROWS, spacing)],
	angles='xy', scale_units='xy', scale=1
)
plt.show()

exit()
#Old version

gravMag = [[0 for j in range(COLS)] for i in range(ROWS)]
grav = [[[0,0] for j in range(COLS)] for i in range(ROWS)]

err = calcDiff(grav, target)
errs = []
for i in range(10000):
	print(i, err)
	errs.append(err)

	randRow = random.randrange(ROWS)
	randCol = random.randrange(COLS)
	randStrength = gravMag[randRow][randCol] + random.gauss(mu=0.0, sigma=0.01)
	if randStrength > 1:
		randStrength = 1
	if randStrength < -1:
		randStrength = -1

	newgrav = deepcopy(grav)
	delta = randStrength - gravMag[randRow][randCol]

	for r in range(ROWS):
		for c in range(COLS):
			dR = r - randRow
			dC = c - randCol
			newgrav[r][c][0] += delta * kernel[dR + ROWS + 1][dC + COLS + 1][0]
			newgrav[r][c][1] += delta * kernel[dR + ROWS + 1][dC + COLS + 1][1]

	newErr = calcDiff(newgrav, target)
	if newErr > err:
		continue
	err = newErr
	grav = newgrav
	gravMag[randRow][randCol] = randStrength

# https://stackoverflow.com/questions/32462881/add-colorbar-to-existing-axis
minimum = min(map(min, gravMag))
maximum = max(map(max, gravMag))

dist = max(abs(minimum), abs(maximum))

fig, axs = plt.subplots(2)
divider = make_axes_locatable(axs[1])
cax = divider.append_axes('right', size='5%', pad=0.05)
axs[0].plot(errs)
im = axs[1].imshow(gravMag, cmap="RdBu", vmin=-dist, vmax=dist)
fig.colorbar(im, cax=cax, orientation='vertical')
plt.show()