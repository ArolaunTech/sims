from copy import deepcopy
import random

"""
Types:
 - 0 = SLTW
 - 1 = CRMC
 - 2 = SPNG/GEL
 - 3 = SNOW/SLTW
"""

rows = 19
cols = 11

def simulate(cell):
	out = 0

	#Fill in powders
	powderCell = deepcopy(cell)
	for r in range(rows):
		for c in range(cols):
			for dc in range(-1, 2):
				nr = r - 1
				nc = c + dc

				if nr == -1 and nc == 0:
					powderCell[r][c] = 3
					break

				if nr < 0 or nr >= rows:
					continue
				if nc < 0 or nc >= cols:
					continue

				if powderCell[nr][nc] == 3: #Powders fall
					powderCell[r][c] = 3
					break

	#Count interactions
	for r in range(rows):
		for c in range(cols):
			if powderCell[r][c] != 3:
				continue
			for dr in range(-1, 2):
				for dc in range(-1, 2):
					nr = r + dr
					nc = c + dc

					if nr < 0 or nr >= rows:
						continue
					if nc < 0 or nc >= cols:
						continue

					if powderCell[nr][nc] == 0:
						out += 1
	return out

curr = [[0 for j in range(cols)] for i in range(rows)]
score = simulate(curr)
for i in range(100000):
	new = deepcopy(curr)
	for j in range(5):
		new[random.randint(0, rows-1)][random.randint(0, cols-1)] = (random.random() < 0.1)

	newScore = simulate(new)
	if newScore > score:
		curr = deepcopy(new)
		score = newScore

	if i % 100 == 0:
		print(i, score)
		for row in curr:
			print("".join([['..', '[]', '<>'][cell] for cell in row]))