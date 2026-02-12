import numpy
import matplotlib.pyplot as plt
import matplotlib.animation as anim

class Plotter:
	def __init__(self, history):
		"""
		'history' is a 3D array where history[t][r][c] is the density at time t, row r, and column c.
		"""
		self.history = history

	def seehist(self):
		fig,ax = plt.subplots()
		ims = []
		for i in range(len(self.history)):
			im = ax.imshow([[self.history[i][r][c] for c in range(len(self.history[0][0]))] for r in range(len(self.history[0]))], animated=True, vmin=0, vmax=1)
			if i == 0:
				ax.imshow([[self.history[i][r][c] for c in range(len(self.history[0][0]))] for r in range(len(self.history[0]))], vmin=0, vmax=1)
			ims.append([im])
		ani = anim.ArtistAnimation(fig, ims, interval = 50, blit = True, repeat_delay = 1000)
		plt.show()

sim = [[0 for i in range(4)] for j in range(4)]
history = [sim]
sim[0][0] = 1
movx = [-1, 0, 1, -1, 0, 1, -1, 0, 1]
movy = [-1, -1, -1, 0, 0, 0, 1, 1, 1]

dt = 0.1

def clamp(x, lo, hi):
	if x > hi:
		return hi
	if x < lo:
		return lo
	return x

for n in range(200):
	new = [[0 for j in range(4)] for i in range(4)]
	for r in range(4):
		for c in range(4):
			for i in range(9):
				nr = r + movx[i]
				nc = c + movy[i]
				nr = clamp(nr, 0, 3)
				nc = clamp(nc, 0, 3)
				new[nr][nc] += dt * sim[r][c] / 9
			new[r][c] += (1 - dt) * sim[r][c]
	sim = new
	history.append(sim)

Plotter(history).seehist()