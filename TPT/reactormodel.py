import matplotlib.pyplot as plt
import numpy as np

def calcReactionTime(initBREL, progresses):
	ticks = 10000
	finalTHRM = progresses*(1-2*abs(initBREL-0.5))
	output = ticks*np.ones(progresses.shape)

	brel = initBREL
	brmt = 1-brel
	thrm = 0

	setInd = 0
	for i in range(ticks):
		THRMgen = 0.001*8*brel*brmt
		brel -= 0.5*THRMgen
		brmt -= 0.5*THRMgen
		thrm += THRMgen

		if thrm > finalTHRM[setInd]:
			output[setInd] = i
			setInd += 1

		#print(brel, brmt, thrm)
	#print(finalTHRM, output)
	return output

#print(calcReactionTime(0.5, np.linspace(0, 1, num=100)))


res=100
initBRELAxis = np.linspace(0,1,num=res)
reactionProgAxis = np.linspace(1,0,num=res)

generationTimes = np.zeros((res, res))
for i in range(res):
	reactionTimes = np.flip(calcReactionTime(initBRELAxis[i], initBRELAxis),axis=0)
	for j in range(res):
		generationTimes[j][i] = reactionTimes[j]

initBREL, reactionProg = np.meshgrid(initBRELAxis, reactionProgAxis)

reactorTHRM = reactionProg*(1-2*np.abs(initBREL-0.5))
reactorBREL = initBREL - 0.5*reactorTHRM
reactorBRMT = 1 - reactorBREL - reactorTHRM
reactorTemp = 3226.85 - 2226.85*np.exp(-20*reactorTHRM)

outputBREL = reactorBREL
outputBRMT = 1 - outputBREL

additionFrac = (1-initBREL)/outputBRMT
use = 1228 + 228*(1-additionFrac)

BRELWaste = (1-additionFrac)*(outputBREL)/initBREL

#plt.pcolor(initBRELAxis, reactionRateAxis, initBREL, cmap='inferno')
outputBREL = np.where(BRELWaste > 0.1, 1, outputBREL)
outputBREL = np.where(reactorTemp < 1860+use, 1, outputBREL)
outputBREL = np.where(generationTimes > 1200, 1, outputBREL)
y_axis = outputBREL
#y_axis = reactorTemp - use
#contour_levels = [0.3,0.5,0.7,0.8,0.9,0.95,0.99]
contour_levels = [1862, 1882, 1901, 1921, 1940]
contour_enabled = True

fig, ax = plt.subplots()

im = ax.imshow(y_axis, cmap='inferno', extent=[0,100,0,100])
if contour_enabled:
	contours = ax.contour(100*initBREL, 100*reactionProg, reactorTemp - use, colors='gray', levels=contour_levels)
	contours2 = ax.contour(100*initBREL, 100*reactionProg, BRELWaste, colors='white', levels=[0.01,0.02,0.03,0.04,0.05])
	contours3 = ax.contour(100*initBREL, 100*reactionProg, generationTimes, colors='red', levels=[100, 200, 300, 400, 500, 600])

ax.set_xlabel('BREL fraction in original mixture (%)')
ax.set_ylabel('THRM generation progress(%)')
ax.set_title('Thermite Reactor Metrics')

cbar = fig.colorbar(im, ax=ax)
#if contour_enabled:
#	cbar.add_lines(contours)
#	cbar.add_lines(contours2)
plt.show()