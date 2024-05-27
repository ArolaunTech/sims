import numpy as np
import matplotlib.pyplot as plt

simWidth = 604
simHeight = 376

simHalfWidthLow = (simWidth//2) - 2
simHalfHeightLow = (simHeight//2) - 2

width = 5
height = 5

def updateSim(simulation):
	pmap = -np.ones((simHeight, simWidth))
	for i in range(len(simulation)):
		pixel = simulation[i]
		pmap[pixel[1]][pixel[2]] = i
	pmap = np.int64(pmap)

	spark = []
	for pixel in simulation:
		if pixel[0] == 1: #PSCN
			if pixel[3] > 0:
				pixel[3] -= 1
		if pixel[0] == 2: #DRAY
			pass
		if pixel[0] == 3: #SPRK
			for rx in range(-2, 3):
				for ry in range(-2, 3):
					x = rx + pixel[1]
					y = ry + pixel[2]

					if x < 0:
						continue
					if y < 0:
						continue
					if x >= simWidth:
						continue
					if y >= simWidth:
						continue
					if pmap[x][y] == -1:
						continue
					if simulation[pmap[x][y]][0] != 1:
						continue
					if simulation[pmap[x][y]][3] != 0:
						continue
					spark.append(pmap[x][y])
			if pixel[3] > 0:
				pixel[3] -= 1
				if pixel[3] == 0:
					pixel[0] = 1
					pixel[3] = 4
	for pixel in spark:
		simulation[pixel][0] = 3
		simulation[pixel][3] = 4
	return simulation

elem_colors = [[0, 0, 0], [128, 80, 80], [255, 170, 34], [255, 255, 128]] #NONE, PSCN, DRAY, SPRK
c = np.array(elem_colors)

while True:
	simulation = [] #Empty field

	for i in range(height): #Add random bomb
		for j in range(width):
			typ = np.random.randint(3)
			if typ == 0:
				continue
			simulation.append([typ, simHalfHeightLow + i, simHalfWidthLow + j, 0])
	simulation = np.array(simulation)

	for i in range(len(simulation)): #For each pixel...
		if simulation[i][0] == 2: #If DRAY, continue
			continue
		simulationTest = np.copy(simulation)
		simulationTest[i][0] = 3 #Spark PSCN
		simulationTest[i][3] = 4
		for i in range(10):
			pmap = np.zeros((simHeight, simWidth))
			for pixel in simulationTest:
				pmap[pixel[1]][pixel[2]] = pixel[0]
			colors = c[np.uint8(pmap)]
			plt.imshow(colors)
			plt.show()

			simulationTest = updateSim(simulationTest)
	print(simulation)