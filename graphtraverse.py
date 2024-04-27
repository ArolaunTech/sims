#Graph traverser
#Includes graph search algs
import math

#Functions
def getInd(node):
	return nodes.index(node)

def setEdge(nodeA,nodeB,cost):
	graph[getInd(nodeA)][getInd(nodeB)] = cost

def setEdgeSymmetric(nodeA, nodeB, cost):
	setEdge(nodeA, nodeB, cost)
	setEdge(nodeB, nodeA, cost)

def minCost(nodeA, nodeB):
	return graph[getInd(nodeA)][getInd(nodeB)]

constd = 315*9.80665*math.log(25/7)
dvstep = 800
dvratio = 1/(1-(1-math.exp(-dvstep/(315*9.80665)))*9/8)
def fuelScore(dv):
	if dv<constd:
		return 14*math.exp(dv/(315*9.80665))
	lastratio = 1/(1-(1-math.exp(-((dv-constd)%dvstep)/(315*9.80665)))*9/8)
	lastratio *= 50
	lastratio *= dvratio ** ((dv-constd)//dvstep)
	return lastratio

#Graph
nodes = [
	'Laythe', 'Laythe Low', 'Laythe High',
	'Vall', 'Vall Low', 'Vall High',
	'Tylo', 'Tylo Low', 'Tylo High',
	'Bop', 'Bop Low', 'Bop High',
	'Pol', 'Pol Low', 'Pol High'
]

graph = [[-1 for node in nodes] for node in nodes]
for i in range(len(nodes)):
	graph[i][i] = 0

#Ascents/Descents
setEdge('Laythe', 'Laythe Low', 2100)
setEdge('Laythe Low', 'Laythe', 0)
setEdgeSymmetric('Vall', 'Vall Low', 900)
setEdgeSymmetric('Tylo', 'Tylo Low', 2400)
setEdgeSymmetric('Bop', 'Bop Low', 250)
setEdgeSymmetric('Pol', 'Pol Low', 150)

#Orbit Changes
setEdge('Laythe Low', 'Laythe High', 596)
setEdge('Laythe High', 'Laythe Low', 0)
setEdgeSymmetric('Tylo Low', 'Tylo High', 798)
setEdgeSymmetric('Vall Low', 'Vall High', 266)
setEdgeSymmetric('Bop Low', 'Bop High', 59)
setEdgeSymmetric('Pol Low', 'Pol High', 41)

#Transfers
setEdgeSymmetric('Pol High', 'Bop High', 197)
setEdgeSymmetric('Pol High', 'Tylo High', 237)
setEdgeSymmetric('Pol High', 'Vall High', 499)
setEdge('Laythe High', 'Pol High', 588)
setEdge('Pol High', 'Laythe High', 385)

setEdge('Laythe High', 'Bop High', 617)
setEdge('Bop High', 'Laythe High', 410)
setEdgeSymmetric('Bop High', 'Vall High', 544)
setEdgeSymmetric('Bop High', 'Tylo High', 360)

setEdge('Laythe High', 'Tylo High', 111)
setEdge('Tylo High', 'Laythe High', 55)
setEdgeSymmetric('Tylo High', 'Vall High', 48)

setEdgeSymmetric('Laythe High', 'Vall High', 33)

#Gravity Assists

setEdgeSymmetric('Laythe High', 'Tylo High', 0) #La-Va-Ty
setEdgeSymmetric('Laythe High', 'Pol High', 179) #La-Va-Ty-Po

setEdgeSymmetric('Laythe High', 'Bop High', 90) #La-Va-Ty-Bo

setEdgeSymmetric('Vall High', 'Bop High', 181) #Va-Ty-Bo
setEdgeSymmetric('Vall High', 'Pol High', 222) #Va-Ty-Ty-Po

print(nodes)
print(graph)

#Propagate graph
for m in range(len(nodes)):
	for i in range(len(nodes)):
		for j in range(len(nodes)):
			if graph[i][j] == -1:
				continue
			for k in range(len(nodes)):
				if graph[j][k] == -1:
					continue
				graph[i][k] = min(graph[i][k], graph[i][j]+graph[j][k])
				if graph[i][k] == -1:
					graph[i][k] = graph[i][j]+graph[j][k]

print(graph)
print('\n')

motherships = ['Laythe Low', 'Laythe High', 'Vall Low', 'Vall High', 'Tylo Low', 'Tylo High', 'Bop Low', 'Bop High', 'Pol Low', 'Pol High']
destinations = ['Laythe', 'Vall', 'Tylo', 'Bop', 'Pol']
departuref = [1, 1, 1, 1, 1]
returnf = [1, 1, 1, 1, 0]

costs = {}

for mothership in motherships:
	costs[mothership] = [minCost(mothership, destinations[i])*departuref[i] + minCost(destinations[i], mothership)*returnf[i] for i in range(len(destinations))]
	print(mothership, costs[mothership])
	for i in range(len(destinations)):
		costs[mothership][i] = fuelScore(costs[mothership][i])
	costs[mothership] = sum(costs[mothership])
print(sorted(zip(costs.keys(),costs.values()), key=lambda x: x[1]))