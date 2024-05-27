import numpy as np
import matplotlib.pyplot as plt

designRows = 9
designCols = 5
zeroCost = 1
misSegmentCost = 1
design = np.zeros((designRows, designCols))

#Set segments on design
for i in range(3):
	for c in [2,3]:
		design[1 + 3*i][c] = i + 5

for r in [2,3]:
	for i in range(2):
		design[r][1 + 3*i] = 1 + 2*i

for r in [5,6]:
	for i in range(2):
		design[r][1 + 3*i] = 2 + 2*i

design[0][0] = 8
design[8][0] = 8

#Functions
def runDesign(design):
	sources = -np.ones((designRows, designCols))
	sources = np.where(np.logical_and(design>0, design<8), 0, sources)

	sensorIndices = np.zeros((designRows, designCols))
	currIndex = 1
	for r in range(len(design)):
		for c in range(len(design[0])):
			if design[r][c] == 9:
				sensorIndices[r][c] = currIndex
				currIndex += 1
	for r in range(len(design)):
		for c in range(len(design[0])):
			if design[r][c] == 9:
				for dr in range(-1,2):
					for dc in range(-1,2):
						nr = (r + dr) % designRows
						nc = (c + dc) % designCols
						if 0 < design[nr][nc] < 8:
							sources[nr][nc] = sensorIndices[r][c]
		if design[r][0] == 9:
			for dr in range(-1,2):
				nr = (r + dr) % designRows
				if 0 < design[nr][-1] < 8:
					sources[nr][-1] = sensorIndices[r][0]
	return sources

def getScore(design, inSources=[]):
	sources = inSources
	if len(inSources) == 0:
		sources = runDesign(design)

	out = 0
	segmentIndices = [[0 for j in range(np.int64(np.max(sources)))] for i in range(7)]
	for r in range(len(sources)):
		for c in range(len(sources[0])):
			if sources[r][c] == 0:
				out += zeroCost
			elif sources[r][c] != -1:
				segmentIndices[np.int64(design[r][c] - 1)][np.int64(sources[r][c] - 1)] += 1
	#print(segmentIndices)

	for j in range(len(segmentIndices[0])):
		maxSource = 0
		sumSource = 0
		for i in range(len(segmentIndices)):
			maxSource = max(segmentIndices[i][j], maxSource)
			sumSource += segmentIndices[i][j]
		sumSource -= maxSource
		sumSource *= misSegmentCost
		out += sumSource
		#print(out)
	return out

def getString(design):
	out = ''
	for r in range(len(design)):
		for c in range(len(design[0])):
			if 0 < design[r][c] < 8:
				out += '\u25A7'
			elif design[r][c] == 9:
				out += '\u25A0'
			else:
				out += '\u25A1'
		out += '\n'
	return out

#Optimization
minScore = 1000000000
optimalDesigns = []
for i in range(10**7):
	design = np.where(design == 9, 0, design)

	for j in range(7):
		addRow = np.random.randint(designRows)
		addCol = np.random.randint(designCols)
		safety = 0
		while design[addRow][addCol] > 0:
			addRow = np.random.randint(designRows)
			addCol = np.random.randint(designCols)
			safety += 1
			if safety > 1000:
				break
		design[addRow][addCol] = 9

	score = getScore(design)

	if score <= minScore:
		if score < minScore:
			optimalDesigns = []
		minScore = score
		alreadyExists = False
		for optimal in optimalDesigns:
			if np.array_equal(optimal, design):
				alreadyExists = True
				break
		if not alreadyExists:
			optimalDesigns.append(design)
			print(score, i, len(optimalDesigns))
			print(getString(design))


#Display
#sources = runDesign(optimalDesigns[0])
#print(getScore(optimalDesigns[0], sources))
#print(len(optimalDesigns))
#fig, ax = plt.subplots(1,2)

#ax[0].imshow(optimalDesigns[0], cmap='inferno')
#ax[1].imshow(sources)
#plt.show()