#Reaction-diffusion system
#Followed tutorial at https://www.algosome.com/articles/reaction-diffusion-gray-scott.html
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import random
import math

class Simulation:
	def __init__(self, r, c, n, reactions, feed, kill, diffusion, state = None):
		#Initializing with state. The state of the grid is a r by c by n array where state[r][c][n] is the concentration of chemical n at (r,c). Chemicals are labeled from 0 to n - 1.
		if state:
			self.state = np.array(state)
		else:
			self.state = np.array([[[0 for k in range(n)] for j in range(c)] for i in range(r)])
		self.r = r
		self.c = c
		self.n = n
		#Reactions list is a 3D array. Each element describes a reaction, where reactions[i][0][n] is the amount of substance n that goes into reaction i, and reactions[i][1][n] is the amount of substance n that comes out of reaction i.
		self.reactions = np.array(reactions)
		#Feed list is a 1D array. feed[n] says how much of substance n you should input into the reaction.
		self.feed = np.array(feed)
		#Kill list is a 1D array. kill[n] says how much of substance n is removed from the reaction.
		self.kill = np.array(kill)
		#Diffusion list is a 1D array. diffusion[n] says how fast substance n diffuses.
		self.diffuse = np.array(diffusion)

	def nabla(self, r, c, state):
		#Calculates diffusion. Given a cell and its neighborhood, it calculates how much substance flows out or into a cell.
		return 0.05 * (state[r-1][c-1] + state[r-1][c+1] + state[r+1][c-1] + state[r+1][c+1]) + 0.2 * (state[r][c-1] + state[r][c+1] + state[r-1][c] + state[r+1][c]) - state[r][c]

	def step(self, state, dt):
		#One step through the simulation.
		new = state
		for i in range(self.r):
			for j in range(self.c):
				if i == 0 or i == self.r - 1 or j == 0 or j == self.c - 1:
					continue
				#For each cell, update it.
				#Diffusion
				new[i][j] += self.diffuse * self.nabla(i, j, state) * dt
				#Reaction
				for reac in self.reactions:
					new[i][j] += np.prod(state[i][j] ** reac[0]) * (reac[1] - reac[0])
				#Feed/Kill
				new[i][j] += self.feed * (1 - state[i][j]) - self.kill * state[i][j]
		return np.clip(new,0,1)

	def run(self, T, dt):
		#Step through the simulation for time T (timestep dt).
		n = math.floor(T / dt)
		history = []
		state = self.state
		for i in range(n):
			history.append(state)
			state = self.step(state, dt)
		return history

class Plotter:
	def __init__(self, history):
		#"history" is a t by r by c by n 4D array where history[t][r][c][n] is the concentration of chemical n at coordinates(r,c) at time t.
		self.history = history

	def show(self, chem):
		fig,ax = plt.subplots()
		ims = []
		for i in range(len(self.history)):
			im = ax.imshow([[self.history[i][r][c][chem] for c in range(len(self.history[0][0]))] for r in range(len(self.history[0]))], animated=True, vmin=0, vmax=1)
			if i == 0:
				ax.imshow([[self.history[i][r][c][chem] for c in range(len(self.history[0][0]))] for r in range(len(self.history[0]))], vmin=0, vmax=1)
			ims.append([im])
		ani = anim.ArtistAnimation(fig, ims, interval = 50, blit = True, repeat_delay = 1000)
		plt.show()

sim = Simulation(100,100,2,[[[1,2],[0,3]]],[0.078,0],[0,0.08],[1,0.5],[[[1,random.random()] for i in range(100)] for j in range(100)])
plot = Plotter(sim.run(500,1))
plot.show(1)