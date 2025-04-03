from utils import *
import random
from copy import deepcopy

g0 = 9.80665

enginesWithFuelTanks = {} #Thrust, fuel/second, fuel/second + tanks, tanks, mass
for engineName in engine_info:
	fuelConsumption = engine_info[engineName][2]/(g0 * engine_info[engineName][0])
	#include fuel tanks
	if engine_info[engineName][1] == "XenonGas":
		fuelConsumption *= 4/3
	if engine_info[engineName][1] == "LiquidFuel":
		fuelConsumption *= 11/10
	if engine_info[engineName][1] == "LFOx":
		fuelConsumption *= 9/8
	if engine_info[engineName][1] == "Monopropellant":
		fuelConsumption *= 8.5/7.5

	enginesWithFuelTanks[engineName] = (
		engine_info[engineName][2], 
		engine_info[engineName][2]/(g0 * engine_info[engineName][0]), 
		fuelConsumption,
		fuelConsumption - engine_info[engineName][2]/(g0 * engine_info[engineName][0]),
		engine_info[engineName][4]
	)

def calculateLeastMassOfStage(payload, minAcceleration, mindV):
	bestConfiguration = {i:0 for i in enginesWithFuelTanks}
	engineConfiguration = deepcopy(bestConfiguration)
	bestMass = 1e9

	for i in range(100000):
		if random.random() < 0.02:
			engineConfiguration = deepcopy(bestConfiguration)
		if random.random() < 0.25:
			engineConfiguration = {i:int(abs(random.gauss(mu=0, sigma=1))) for i in enginesWithFuelTanks}
		else:
			changeEngine = random.choice(list(engineConfiguration.keys()))
			if random.random() < 0.5:
				engineConfiguration[changeEngine] += 1
			elif engineConfiguration[changeEngine] > 0:
				engineConfiguration[changeEngine] -= 1

		totalThrust = 0
		totalFuelConsumption = 0
		totalTankAddition = 0
		totalDryMass = payload
		for engine in engineConfiguration:
			totalThrust += enginesWithFuelTanks[engine][0] * engineConfiguration[engine]
			totalFuelConsumption += enginesWithFuelTanks[engine][1] * engineConfiguration[engine]
			totalTankAddition += enginesWithFuelTanks[engine][3] * engineConfiguration[engine]
			totalDryMass += enginesWithFuelTanks[engine][4] * engineConfiguration[engine]

		if totalThrust == 0:
			continue

		veff = totalThrust/totalFuelConsumption
		massratio = math.exp(mindV/veff)

		fuelmult = (totalDryMass * totalFuelConsumption/(totalFuelConsumption + totalTankAddition - massratio * totalTankAddition) - totalDryMass)/totalTankAddition
		wetmass = totalDryMass + fuelmult * (totalTankAddition + totalFuelConsumption)
		acceleration = totalThrust/wetmass

		if acceleration < minAcceleration:
			continue

		if wetmass < bestMass:
			bestConfiguration = deepcopy(engineConfiguration)
			bestMass = wetmass

	#print(bestConfiguration)
	for engine in bestConfiguration:
		if bestConfiguration[engine] > 0:
			print(f"{bestConfiguration[engine]}x {engine}")
	return bestMass

print(calculateLeastMassOfStage(1000, 9, 740))