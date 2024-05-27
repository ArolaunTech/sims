#Fuel optimizer
import math
import random

#Constants
g = 9.80665

#Functions
def dv(isp, w, d):
	return isp*g*math.log(w/d)

def calcscore(fuels, required):
	s = sum(fuels)
	actfuel = [fuels[i]*required/s for i in range(len(fuels))]
	out = 0
	for i in range(1, len(fuels)+1):
		wet = p+pt*i+sum(actfuel[:i])+w*i/len(fuels)
		dry = p+pt*i+sum(actfuel[:i-1])+w*i/len(fuels)+actfuel[i-1]/9

		out += dv(cisp, wet, dry)
	return out

def closestwhole(distribution):
	newdist = [math.floor(distribution[i]) for i in range(len(distribution))]
	frac = [distribution[i] % 1 for i in range(len(distribution))]
	sumfrac = round(sum(frac))
	det = sorted([i for i in range(len(distribution))], key=lambda x:frac[x], reverse=True)
	for i in range(sumfrac):
		newdist[det[i]] += 1
	return newdist


#Craft constants
mass = 50 #Total mass in tons
p = (4+0.03)+(0.16+0.1+0.01+0.08+0.05+0.09) #Payload mass
p -= 0.025*2 #Correct for decoupler
cisp = 295

pt = 2*((0.075+0.03)+0.025)
w = 1.2+0.08*2

for T in range(4, 10):

	#Calculations
	fuelm = 50-p-w-pt*T

	dist = [1/T for i in range(T)]

	for i in range(1000):
		i1 = random.randint(0, T-1)
		i2 = random.randint(0, T-2)

		if i2 >= i1:
			i2 += 1

		c = random.random()
		c *= 0.05*min(dist[i1],dist[i2])

		score = calcscore(dist, fuelm)

		dist[i1] -= c
		dist[i2] += c

		if score >= calcscore(dist, fuelm) or dist[i1] < 4/fuelm:
			dist[i1] += c
			dist[i2] -= c
	print(score, closestwhole([dist[i] * fuelm for i in range(T)]), calcscore(closestwhole([dist[i] * fuelm for i in range(T)]), fuelm))
	print(fuelm, T, calcscore([1,1,1],fuelm))
