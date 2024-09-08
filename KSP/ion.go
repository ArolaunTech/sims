package main

import (
	"fmt"
	"math/rand"
	"math"
)

const g float64 = 9.81
const kerbinRadius float64 = 6.0e+5
const kerbinStdgp float64 = g * kerbinRadius * kerbinRadius
const kerbinRotp float64 = 21549.425
const kerbinAtmoAlt float64 = 7.0e+4

const maxKerbinCircAlt float64 = 9.0e+6
const minKerbinCircAlt float64 = 7.0e+4

const spaceKerbinAlt float64 = 7.0e+4
const highKerbinAlt float64 = 1.8e+4

const mohoRadius float64 = 2.5e+5
const mohoStdgp float64 = 0.275 * g * mohoRadius * mohoRadius
const mohoRotp float64 = 1.21e+6

const mohoOrbAlt float64 = 1.0e+4
const mohoLandAlt float64 = 5.0e+3

const mohoRV float64 = 2.0 * math.Pi * (mohoRadius + mohoLandAlt)/mohoRotp
var mohoGAccel float64 = mohoStdgp/math.Pow(mohoRadius + mohoLandAlt, 2.0)
var mohoLandOrbVel float64 = math.Sqrt(mohoStdgp/(mohoRadius + mohoLandAlt))

const ionEC float64 = 600.0/(7.0 * g)
const ionXe float64 = 1.0/(2100.0 * g)

const rapierIspClosedCycle float64 = 305.0
const rapierClosedCycleThrust float64 = 180.0
const rapierMaxVel float64 = 1680.0
const rapierMass float64 = 2

func calcSingleEngineMax(ec float64) (float64, float64, float64) {
	//Returns single ion engine thrust (kN), total EC consumption (u/s), and xenon consumption (tons/s).

	if ec <= 0.0 {
		return 0.0, 0.0, 0.0
	}
	if ec >= ionEC {
		return 2.0, ionEC, 2.0/(4200.0 * g)
	}
	if ec <= 0.1 * ionEC {
		return 11.0 * ec/(ionEC), ec, 2.0 * ec/(4200.0 * g * ionEC * 0.1)
	}

	xenonConsump := 2.0 / (4200.0 * g)
	thrustConsump := 1.0 + ec/ionEC
	return thrustConsump, ec, xenonConsump
}

func calcIonStarvation(ionThrust float64, throttle float64, currEC float64, solar float64) (float64, float64, float64) {
	//Calculates ion engine thrust, total EC consumption, and xenon consumption.
	//Returns:
	// - Thrust (float64, kN)
	// - EC Consumption (float64, u/s)
	// - Xenon Consumption (float64, tons/s)

	numEngines := ionThrust/2.0

	demandXe := ionThrust * throttle/(4200.0 * g)
	demandEC := ionEC * numEngines * throttle

	if (currEC > 0) || (solar >= demandEC) { //No ion starvation
		return ionThrust*throttle, demandEC, demandXe
	}
	if solar <= 0 {
		return 0.0, 0.0, 0.0
	}

	//Here it is guaranteed that the ion engines are not getting enough EC to run
	desiredStarve := 0.1 * ionEC * throttle
	if solar <= desiredStarve {
		return 1.1 * solar * throttle/desiredStarve, solar, ionXe * throttle * solar/desiredStarve
	}
	desiredFull := demandEC - ionEC * throttle
	desiredFullThrust := desiredFull * 2.0/ionEC
	desiredFullXe := desiredFullThrust/(4200.0 * g)

	if solar >= desiredFull {
		lastThrust, _, lastXe := calcSingleEngineMax((solar - desiredFull)/throttle)
		return desiredFullThrust + lastThrust*throttle, solar, desiredFullXe + lastXe*throttle
	}
	return 1.1*throttle + 2.0*(solar-desiredStarve)/ionEC, solar, ionXe*(throttle+(solar-desiredStarve)/ionEC)
}

func neededLF(mass float64) float64 {
	lf := mass * 0.065
	if lf > 0.5 {
		return 0.5 + math.Ceil((lf - 0.5) * 40) * 0.025
	}
	return math.Ceil(lf * 20.0) * 0.05
}

func minXeCapacity(xenon float64) float64 {
	if xenon <= 0.0 {
		return 0.0
	}
	//if xenon >= 0.642 {
	//	return 0.642 + minXeCapacity(xenon - 0.642)
	//}
	out := 0.0405 + minXeCapacity(xenon - 0.0405)
	out = math.Min(out, 0.072 + minXeCapacity(xenon - 0.072))
	out = math.Min(out, 0.57 + minXeCapacity(xenon - 0.57))
	return out
}

func minLFOXCapacity(lfox float64) float64 {
	if lfox <= 0.0 {
		return 0.0
	}
	out := 0.11 + minLFOXCapacity(lfox - 0.11)
	out = math.Min(out, 0.27 + minLFOXCapacity(lfox - 0.27))
	out = math.Min(out, 0.2 + minLFOXCapacity(lfox - 0.2))
	out = math.Min(out, 0.3 + minLFOXCapacity(lfox - 0.3))
	out = math.Min(out, 0.5 + minLFOXCapacity(lfox - 0.5))
	return out
}

func simulateMohoLandingAscent(mass float64, battery float64, solar float64, ionThrust float64, chemThrust float64) (float64, float64, float64, float64) {
	//Returns mass, punishMass, xenon, LFOX
	//nextCodeIndex := 0
	currCodeIonThrust := 1.0
	currCodeChemThrust := 0.0
	//currCodeTime := 0.0
	currMass := mass
	currEC := battery
	totalLFOXUsed := 0.0
	totalXenonUsed := 0.0

	maxChemLFOXConsumption := chemThrust/(rapierIspClosedCycle * g)

	//var ionPush, ionXeConsump, ionECConsump float64
	var chemPush, chemLFOXConsump, totalPush float64

	var accel, vaccel, haccel float64

	lowMaxSpeed := math.Sqrt(mohoStdgp*(2.0/(mohoRadius + mohoLandAlt) - 2.0/(2.0 * mohoRadius + mohoLandAlt + mohoOrbAlt)))
	currVel := lowMaxSpeed

	ascentValid := false
	//codeTime := 0.0

	timestep := 0.5
	for time := 0.0; time < 3000.0; time += timestep {
		if math.Mod(time, 1.0) == 0.0 {
			fmt.Printf("T = %fs: Velocity of %f m/s, Currently %fkg, XE Used: %fkg, LF+OX Used: %fkg\n", time, currVel, 1000.0*currMass, 1000.0*totalXenonUsed, 1000.0*totalLFOXUsed)
		}
		//Read code
		//if nextCodeIndex < len(code) {
		//	if time >= currCodeTime + code[nextCodeIndex] {
		//		currCodeTime += code[nextCodeIndex]
		//		currCodeChemThrust = code[nextCodeIndex + 2]
		//		nextCodeIndex += 3
		//	}
		//}

		//Calculate thrust
		vaccel = mohoGAccel * math.Abs(1.0 - math.Pow(currVel/mohoLandOrbVel,2.0))

		ionPush, iEC, ionXeConsump := calcIonStarvation(ionThrust, currCodeIonThrust, currEC, solar)
		netEC := solar - iEC

		currCodeChemThrust = (1.125 * maxChemLFOXConsumption * ionPush * ionPush - 1.125 * maxChemLFOXConsumption * currMass * currMass * vaccel * vaccel - 4.0/3.0 * ionXeConsump * ionPush * chemThrust)/(4.0/3.0 * ionXeConsump * chemThrust * chemThrust - 1.125 * maxChemLFOXConsumption * chemThrust * ionPush)
		currCodeChemThrust = math.Max(currCodeChemThrust, 0.0)
		currCodeChemThrust = math.Min(currCodeChemThrust, 1.0)

		chemPush = currCodeChemThrust * chemThrust
		chemLFOXConsump = currCodeChemThrust * maxChemLFOXConsumption

		totalPush = chemPush + ionPush

		//Calculate acceleration
		accel = totalPush / currMass
		if accel < vaccel {
			return -1.0, -1.0, -1.0, -1.0
		}
		haccel = math.Sqrt(accel * accel - vaccel * vaccel)

		//Updates
		currMass -= chemLFOXConsump * timestep
		currMass -= ionXeConsump * timestep
		currEC += timestep * netEC
		currEC = math.Max(currEC, 0.0)
		currEC = math.Min(currEC, battery)

		totalLFOXUsed += chemLFOXConsump * timestep
		totalXenonUsed += ionXeConsump * timestep

		currVel -= haccel * timestep
		if currVel <= mohoRV {
			ascentValid = true
			//codeTime = time
			break
		}

		if currMass <= 0.0 {
			return -1.0, -1.0, -1.0, -1.0
		}

		//fmt.Print(currVel, ionPush, chemPush, totalPush, currMass, currEC, "\n")
	}

	if !ascentValid {
		return -1.0, -1.0, -1.0, -1.0
	}

	currVel = mohoRV
	currEC = battery

	for time := 0.0; time < 3000.0; time += timestep {
		if math.Mod(time, 1.0) == 0.0 {
			//fmt.Print(time, currVel, totalXenonUsed, totalLFOXUsed, currMass, "\n")
			fmt.Printf("T = %fs: Velocity of %f m/s, Currently %fkg, XE Used: %fkg, LF+OX Used: %fkg\n", time, currVel, 1000.0*currMass, 1000.0*totalXenonUsed, 1000.0*totalLFOXUsed)
		}
		//Read code
		//if nextCodeIndex < len(code) {
		//	if time + codeTime >= currCodeTime + code[nextCodeIndex] {
		//		currCodeTime += code[nextCodeIndex]
		//		currCodeChemThrust = code[nextCodeIndex + 2]
		//		nextCodeIndex += 3
		//	}
		//}

		//Calculate thrust
		vaccel = mohoGAccel * math.Abs(1.0 - math.Pow(currVel/mohoLandOrbVel,2.0))

		ionPush, iEC, ionXeConsump := calcIonStarvation(ionThrust, currCodeIonThrust, currEC, solar)
		netEC := solar - iEC

		currCodeChemThrust = (1.125 * maxChemLFOXConsumption * ionPush * ionPush - 1.125 * maxChemLFOXConsumption * currMass * currMass * vaccel * vaccel - 4.0/3.0 * ionXeConsump * ionPush * chemThrust)/(4.0/3.0 * ionXeConsump * chemThrust * chemThrust - 1.125 * maxChemLFOXConsumption * chemThrust * ionPush)
		currCodeChemThrust = math.Max(currCodeChemThrust, 0.0)
		currCodeChemThrust = math.Min(currCodeChemThrust, 1.0)

		chemPush = currCodeChemThrust * chemThrust
		chemLFOXConsump = chemPush/(rapierIspClosedCycle * g)

		totalPush = chemPush + ionPush

		//Calculate acceleration
		accel = totalPush / currMass
		if accel < vaccel {
			return -1.0, -1.0, -1.0, -1.0
		}
		haccel = math.Sqrt(accel * accel - vaccel * vaccel)

		//Updates
		currMass -= chemLFOXConsump * timestep
		currMass -= ionXeConsump * timestep
		currEC += timestep * netEC
		currEC = math.Max(currEC, 0.0)
		currEC = math.Min(currEC, battery)

		totalLFOXUsed += chemLFOXConsump * timestep
		totalXenonUsed += ionXeConsump * timestep

		currVel += haccel * timestep
		if currVel >= lowMaxSpeed {
			currPMass := currMass - (totalXenonUsed/3.0) - (totalLFOXUsed/8.0)
			return currMass, currPMass, totalXenonUsed, totalLFOXUsed
		}

		if currMass <= 0.0 {
			return -1.0, -1.0, -1.0, -1.0
		}

		//fmt.Print(currVel, ionPush, chemPush, totalPush, currMass, currEC, "\n")
	}

	return -1.0, -1.0, -1.0, -1.0
}

func simulateIonCircularization(apoapsis float64, periapsis float64, code []float64, mass float64, battery float64, solar float64, ionThrust float64, chemThrust float64, fcConsump float64) (float64, float64, float64, float64, float64, float64) {
	startingAlt := 670000.0
	sma := 0.5*(apoapsis + periapsis)
	startingVel := math.Sqrt(kerbinStdgp * (2.0/startingAlt - 1.0/sma))

	eccentricity := apoapsis/sma - 1.0
	horizontalVel := math.Sqrt((1.0 - math.Pow(eccentricity, 2.0))*sma*kerbinStdgp)/startingAlt
	verticalVel := math.Sqrt(math.Pow(startingVel, 2.0) - math.Pow(horizontalVel, 2.0))

	//fmt.Printf("%s%f%s\n%s%f%s\n", "Hvel: ", horizontalVel, "m/s", "Vvel: ", verticalVel, "m/s")

	x := 0.0
	y := startingAlt
	vx := horizontalVel
	vy := verticalVel
	currMass := mass
	currEC := battery

	nextCodeIndex := 0
	currCodeAngle := 0.0
	currCodeIonThrust := 1.0
	currCodeChemThrust := 0.0
	totalXenonUsed := 0.0
	totalLFOXUsed := 0.0

	currCodeTime := 0.0

	timestep := 0.5

	var distance, gx, gy, gForce float64

	//var ionPush, ionXeConsump, ionECConsump float64
	var chemPush, chemLFOXConsump, totalPush float64

	var zX, zY float64
	var craftX, craftY float64
	var aC, aS float64

	var tX, tY float64

	var currSMA, currEccen, currPeri, currApo float64
	var periVel, massRewarded, massPunished float64

	//var zDist float64

	for time := 0.0; time < 3000.0; time += timestep {
		//Calculate gravity
		distance = math.Hypot(x, y)
		gForce = kerbinStdgp/math.Pow(distance, 3.0)
		gx = -x * gForce
		gy = -y * gForce

		//Read code
		if nextCodeIndex < len(code) {
			if time > code[nextCodeIndex] + currCodeTime {
				currCodeTime += code[nextCodeIndex]
				currCodeAngle = code[nextCodeIndex + 1]
				currCodeIonThrust = code[nextCodeIndex + 2]
				currCodeChemThrust = code[nextCodeIndex + 3]
				nextCodeIndex += 4

				aC, aS = math.Cos(currCodeAngle), math.Sin(currCodeAngle)
			}
		}

		//Calculate thrust
		ionPush, iEC, ionXeConsump := calcIonStarvation(ionThrust, currCodeIonThrust, currEC, solar)
		netEC := solar - iEC
		//fmt.Print(ionPush, "\n")

		chemPush = chemThrust * currCodeChemThrust
		chemLFOXConsump = chemPush/(rapierIspClosedCycle*g)
		totalPush = chemPush + ionPush

		//Calculate craft direction
		zX = y/distance
		zY = -x/distance
		//zX = 1
		//zY = 0

		craftX = zX * aC - zY * aS
		craftY = zY * aC + zX * aS

		//Calculate thrust components
		tX, tY = craftX * totalPush, craftY * totalPush
		//fmt.Print(tX/currMass, tY/currMass, gx, gy, "\n")

		//Update velocities
		vx += gx * timestep
		vy += gy * timestep

		vx += tX * timestep/currMass
		vy += tY * timestep/currMass

		//Update positions
		x += vx * timestep
		y += vy * timestep

		//Other updates
		currMass -= chemLFOXConsump * timestep
		currMass -= fcConsump * timestep
		totalLFOXUsed += chemLFOXConsump * timestep
		totalLFOXUsed += fcConsump * timestep

		currMass -= ionXeConsump * timestep
		totalXenonUsed += ionXeConsump * timestep

		currEC += netEC * timestep
		currEC = math.Min(currEC, battery)
		currEC = math.Max(currEC, 0.0)

		//Calculate Orbital Parameters
		currSMA = 1.0/(2.0/distance - (math.Pow(vx, 2.0) + math.Pow(vy, 2.0))/kerbinStdgp)
		currEccen = math.Sqrt(1.0 - math.Pow(vx * y - vy * x, 2.0)/(currSMA * kerbinStdgp))
		currPeri = currSMA * (1.0 - currEccen)
		currApo = currSMA * (1.0 + currEccen)
		//fmt.Print(currPeri, currApo, "\n")

		if distance < kerbinRadius + minKerbinCircAlt {
			return -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
		}
		if distance > kerbinRadius + maxKerbinCircAlt {
			return -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
		}
		if currApo > 9.0e+6 {
			return -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
		}
		if currMass < 0.0 {
			return -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
		}

		if currPeri >= kerbinRadius + spaceKerbinAlt {
			periVel = math.Sqrt(kerbinStdgp * (2.0/currPeri - 1.0/currSMA)) - math.Sqrt(kerbinStdgp/currPeri)
			massRewarded = currMass * math.Exp(periVel/(4200.0 * g))
			massPunished = massRewarded - (totalXenonUsed + currMass - massRewarded)/3.0 - totalLFOXUsed/8.0
			return currMass, massRewarded, massPunished, totalXenonUsed, totalLFOXUsed, currApo
		}
		//fmt.Print(time, distance, x, y, vx, vy, currMass, currEC, currCodeAngle, currCodeIonThrust, currCodeChemThrust, totalPush, "\n")
	}

	return -1.0, -1.0, -1.0, -1.0, -1.0, -1.0
}

func calcOptimalKerbinAscent(mass float64, solar float64, battery float64, ionThrust float64, chemThrust float64, fcConsump float64) (float64, float64, float64, float64, float64, float64) {
	bestMass, bestRewardMass, bestPunishMass := -1.0, -1.0, -1.0
	bestXenon, bestLFOX, bestLF := -1.0, -1.0, -1.0
	bestApoapsis := -1.0
	numImprovements := 0

	bestCode := []float64{}
	bestSApo := -1.0
	bestSPeri := -1.0

	//Calculate the LF-only part of flight
	bestLF = neededLF(mass)
	highAltMass := mass - bestLF

	var randMass, randRewardMass, randPunishMass, randXe, randLFOX, randApo float64
	var highVel, highHVel, highVVel, highEccen, highSMA, closedCycleMass float64

	for _ = range 50 {
		bestRandomMass, bestRandomRewardMass, bestRandomPunishMass, bestRandomXe, bestRandomLFOX, bestRandomApo, bestRandomApoRes, bestRandomPeri, bestRandomCode := -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, []float64{}
		var randomApo, randomPeri float64
		for iterRand := range 2000 {
			//Generate random apoapsis, periapsis, and codes
			randomApo = kerbinRadius + kerbinAtmoAlt + 45000.0 + 20000.0 * rand.Float64()
			randomPeri = (kerbinRadius + highKerbinAlt) * (rand.Float64() * 0.3 + 0.2)
			randomCode := []float64{}
			for _ = range 1 + rand.Intn(4) {
				randomCode = append(randomCode, 300.0 * rand.Float64(), math.Pi * (rand.Float64() - 0.5), 1.0, rand.Float64())
			}

			if iterRand < 750 {
				randomCode = []float64{0.0, 0.3 + 0.1*rand.Float64(), 1.0, 0.0}
			}

			highVel = math.Sqrt(kerbinStdgp * (2/(kerbinRadius + highKerbinAlt) - 2/(randomApo + randomPeri)))
			highSMA = 0.5 * (randomApo + randomPeri)
			highEccen = randomApo/highSMA - 1.0
			highHVel = math.Sqrt((1 - math.Pow(highEccen, 2.0)) * highSMA * kerbinStdgp)/(kerbinRadius + highKerbinAlt)
			highVVel = math.Sqrt(highVel * highVel - highHVel * highHVel)

			//The planet rotates beneath us
			//fmt.Print(highVel, highSMA, highEccen, highHVel, highVVel, "\n")
			highHVel -= 2 * math.Pi * (kerbinRadius + highKerbinAlt)/kerbinRotp

			//RAPIER burning
			highVel = math.Hypot(highHVel, highVVel)
			highVel -= rapierMaxVel
			if highVel < 0.0 {
				highVel = 0.0
			}

			//Calculate mass at the end of closed cycle burning
			closedCycleMass = highAltMass/math.Exp(highVel/(rapierIspClosedCycle * g))

			//fmt.Print(highVel, highHVel, highVVel, closedCycleMass, "\n\n")

			randMass, randRewardMass, randPunishMass, randXe, randLFOX, randApo = simulateIonCircularization(randomApo, randomPeri, randomCode, closedCycleMass, battery, solar, ionThrust, chemThrust, fcConsump)
			if randMass > -1.0 {
				//fmt.Print(highVel, highHVel, highVVel, closedCycleMass, "\n\n")
				//fmt.Print(randMass, randRewardMass, randPunishMass, randXe, randLFOX, randApo, "\n\n")
				//fmt.Print(randomApo, randomPeri, randomCode, "\n\n")

				if randPunishMass > bestRandomPunishMass {
					bestRandomMass = randMass
					bestRandomRewardMass = randRewardMass
					bestRandomPunishMass = randPunishMass
					bestRandomXe = randXe
					bestRandomLFOX = randLFOX+highAltMass - closedCycleMass
					bestRandomApoRes = randApo
					bestRandomApo = randomApo
					bestRandomPeri = randomPeri
					bestRandomCode = randomCode
				}
			}
		}
		if bestRandomMass == -1.0 {
			continue
		}

		//Iteratively improve the code to make it better over time
		if (bestRandomPunishMass > bestPunishMass) {fmt.Printf("New Kerbin ascent found:\n - Circularization code: %v\n - Starting periapsis: %fm\n - Starting apoapsis: %fm\n - Final mass: %ft\n - Final apoapsis: %fm\n\n", bestRandomCode, bestRandomPeri - kerbinRadius, bestRandomApo - kerbinRadius, bestRandomMass, bestRandomApoRes - kerbinRadius)}
		
		oldCode := bestRandomCode
		oldApo := bestRandomApo
		oldPeri := bestRandomPeri

		oldMass := bestRandomMass
		oldRMass := bestRandomRewardMass
		oldPMass := bestRandomPunishMass
		oldXe := bestRandomXe
		oldLFOX := bestRandomLFOX
		oldFApo := bestRandomApoRes

		apoVel, periVel, codeVel := 0.0, 0.0, []float64{}
		for _ = range len(oldCode) {
			codeVel = append(codeVel, 0.0)
		}

		newCode := oldCode
		newApo, newPeri := oldApo, oldPeri
		for iter := range 10000 {
			apoVel *= 0.9
			periVel *= 0.9

			apoVel += 1.0 * rand.NormFloat64()
			periVel += 1.0 * rand.NormFloat64()

			for i := 0; i < len(oldCode); i += 4 {
				//fmt.Print(i)
				codeVel[i] *= 0.9
				codeVel[i + 1] *= 0.9
				codeVel[i + 2] *= 0.9
				codeVel[i + 3] *= 0.9

				codeVel[i] += 10000.0 * rand.NormFloat64() + 0.1 * (newCode[i] - oldCode[i])
				codeVel[i + 1] += 0.5 * rand.NormFloat64() + 0.1 * (newCode[i + 1] - oldCode[i + 1])
				codeVel[i + 2] += 0.5 * rand.NormFloat64() + 0.1 * (newCode[i + 2] - oldCode[i + 2])
				codeVel[i + 3] += 0.5 * rand.NormFloat64() + 0.1 * (newCode[i + 3] - oldCode[i + 3])
			}

			for i := 0; i < len(oldCode); i += 4 {
				newCode[i] += codeVel[i]
				newCode[i + 1] += codeVel[i + 1]
				newCode[i + 2] += codeVel[i + 2]
				newCode[i + 3] += codeVel[i + 3]

				newCode[i] = math.Max(newCode[i], 0.0)

				newCode[i + 1] = math.Min(newCode[i + 1], math.Pi * 0.5)
				newCode[i + 1] = math.Max(newCode[i + 1], -math.Pi * 0.5)

				newCode[i + 2] = math.Min(newCode[i + 2], 1.0)
				newCode[i + 2] = math.Max(newCode[i + 2], 0.0)

				newCode[i + 3] = math.Min(newCode[i + 3], 1.0)
				newCode[i + 3] = math.Max(newCode[i + 3], 0.0)
			}
			newApo += apoVel
			newPeri += periVel

			newApo = math.Min(newApo, kerbinRadius + maxKerbinCircAlt)
			newApo = math.Max(newApo, kerbinRadius + kerbinAtmoAlt)

			newPeri = math.Min(newPeri, kerbinRadius + highKerbinAlt)
			newPeri = math.Max(newPeri, 0.0)

			//Simulation
			highVel = math.Sqrt(kerbinStdgp * (2/(kerbinRadius + highKerbinAlt) - 2/(newApo + newPeri)))
			highSMA = 0.5 * (newApo + newPeri)
			highEccen = newApo/highSMA - 1.0
			highHVel = math.Sqrt((1 - math.Pow(highEccen, 2.0)) * highSMA * kerbinStdgp)/(kerbinRadius + highKerbinAlt)
			highVVel = math.Sqrt(highVel * highVel - highHVel * highHVel)

			//The planet rotates beneath us
			//fmt.Print(highVel, highSMA, highEccen, highHVel, highVVel, "\n")
			highHVel -= 2 * math.Pi * (kerbinRadius + highKerbinAlt)/kerbinRotp

			//RAPIER burning
			highVel = math.Hypot(highHVel, highVVel)
			highVel -= rapierMaxVel
			if highVel < 0.0 {
				highVel = 0.0
			}

			//Calculate mass at the end of closed cycle burning
			closedCycleMass = highAltMass/math.Exp(highVel/(rapierIspClosedCycle * g))

			newMass, newRMAss, newPMass, newXe, newLFOX, newFApo := simulateIonCircularization(newApo, newPeri, newCode, closedCycleMass, battery, solar, ionThrust, chemThrust, fcConsump)
			newLFOX += highAltMass - closedCycleMass

			if newPMass > -1.0 {
				if newPMass > oldPMass {
					oldCode, oldApo, oldPeri, oldPMass, oldMass, oldRMass, oldXe, oldLFOX, oldFApo = newCode, newApo, newPeri, newPMass, newMass, newRMAss, newXe, newLFOX, newFApo
					numImprovements += 1
					if newPMass > bestPunishMass {
						fmt.Printf("Kerbin ascent found on iteration #%d:\n - Circularization code: %v\n - Starting periapsis: %fm\n - Starting apoapsis: %fm\n - Final mass: %ft\n - Final apoapsis: %fm\n\n", iter, newCode, newPeri - kerbinRadius, newApo - kerbinRadius, newMass, newFApo - kerbinRadius)
						//fmt.Print(newMass, iter, newApo, newPeri, newFApo, "\n", newCode, "\n\n")
					}
				} else {
					apoVel *= 0.9
					periVel *= 0.9
					for i := range codeVel {
						codeVel[i] *= 0.9
					}
				}
			} else {
				apoVel *= 0.7
				periVel *= 0.7
				for i := range codeVel {
					codeVel[i] *= 0.7
				}
			}

			if iter % 200 == 0 {
				newCode, newApo, newPeri = oldCode, oldApo, oldPeri
				apoVel = 0.0
				periVel = 0.0
				for i := range codeVel {
					codeVel[i] = 0.0
				}
			}
		}

		if oldPMass > bestPunishMass {
			bestMass = oldMass
			bestRewardMass = oldRMass
			bestPunishMass = oldPMass
			bestXenon = oldXe
			bestLFOX = oldLFOX
			bestApoapsis = oldFApo
			bestCode = oldCode
			bestSApo = oldApo
			bestSPeri = oldPeri
			//fmt.Print(oldPMass, "\n")
		}
	}

	fmt.Print(bestCode, bestSApo, bestSPeri, "\n\n")

	//bestMass += 1.0
	//bestRewardMass += 1.0
	//bestPunishMass += 1.0
	//bestXenon += 1.0
	//bestLFOX += 1.0
	//bestApoapsis += 1.0

	//fmt.Print(numImprovements, "\n")

	return bestMass, bestRewardMass, bestXenon, bestLFOX, bestLF, bestApoapsis //mass, massReward, xenon, lfox, lf, apoapsis
}

func main() {
	//t0, t1, t2 := calcIonStarvation(14.0, 1.0, 0.0, 11.48)
	//fmt.Print(t0, t1, t2, "\n\n")

	//fmt.Print(neededLF(6.87), "\n\n")

	//mass, rewardMass, punishMass, xenon, lfox, apo := simulateIonCircularization(697628, 274915, []float64{0.0, 0.217, 1.0, 0.0}, 5.5, 1000.0, 4.92, 14.0, 180.0)
	//fmt.Print(mass, rewardMass, punishMass, xenon, lfox, apo, "\n\n")
	t, ec, xe := calcIonStarvation(14.0, 1.0, 0.0, 10.0)
	fmt.Print(t, ec, xe, "\n")
	t, ec, xe = calcIonStarvation(14.0, 0.5, 0.0, 10.0)
	fmt.Print(t, ec, xe, "\n\n")

	numIons := 7.0
	numPanels := 0.0
	numBigPanels := 1.0
	batteryCapacity := 100.0
	startingMass := 6.61
	fcArrays := 0.0

	fmt.Printf("Starting Mass: %ft\n\n", startingMass)

	mass, rewardMass, xenon, lfox, lf, apo := calcOptimalKerbinAscent(startingMass, 1.64 * numPanels + 24.4 * numBigPanels + 18.0 * fcArrays, batteryCapacity, 2.0 * numIons, rapierClosedCycleThrust, 0.000225 * fcArrays)
	fmt.Printf("Ascent Stats:\n - Final Mass: %ft\n - Xenon Used: %ft\n - LF+OX mixture used: %ft\n - LF used: %ft\n - Final Apoapsis: %fm\n\n",mass, xenon, lfox, lf, apo - kerbinRadius)
	mohoMass := rewardMass/math.Exp(1500.0/(4200.0*g))
	fmt.Printf("Mass at Moho:\n%ft\n\n", mohoMass)

	//orbMass, _, orbXe, orbLFOX := calcOptimalMohoLanding(mohoMass, batteryCapacity, 16.4 * numPanels + 244.0 * numBigPanels, 2.0 * numIons, rapierClosedCycleThrust)
	orbMass, _, orbXe, orbLFOX := simulateMohoLandingAscent(mohoMass, batteryCapacity, 16.4 * numPanels + 244.0 * numBigPanels, 2.0 * numIons, rapierClosedCycleThrust)
	totalLFOX := orbLFOX + lfox
	totalXe := orbXe + xenon + mass - mohoMass

	dryMass := orbMass/math.Exp(600.0/(4200*g))
	totalXe += orbMass - dryMass

	fmt.Printf("Final Mass: %ft\n\n", dryMass)

	fmt.Printf("Total xenon used: %ft (%f units)\nTotal LF+OX mixture used: %ft (%f units)\nTotal LF used: %ft (%f units)\n\n", totalXe, 10000.0 * totalXe, totalLFOX, 200.0 * totalLFOX, lf, 200.0 * lf)
	//fmt.Print(totalXe, totalLFOX, lf, "\n")
	actualOX := totalLFOX * 0.55
	actualLF := totalLFOX * 0.45 + lf

	fmt.Printf("Real LF used: %ft (%f units)\nReal OX used: %ft (%f units)\n\n", actualLF, 200.0 * actualLF, actualOX, 200.0 * actualOX)
	//fmt.Print(10000.0 * totalXe, 200.0 * actualOX, 200.0 * actualLF, "\n")

	minimumDryMass := 0.0
	minimumDryMass += rapierMass //RAPIER
	minimumDryMass += 0.25 * numIons //Ion engines
	minimumDryMass += 0.2 //Fairing
	minimumDryMass += 0.1 //Wing
	minimumDryMass += 0.025 * numPanels //Solar Panels
	minimumDryMass += 0.3 * numBigPanels //Gigantors
	minimumDryMass += 0.14 //Kerbal + seat
	minimumDryMass += 0.05 //Reaction wheel
	minimumDryMass += batteryCapacity * 0.00005 //Batteries
	minimumDryMass += 0.02 //Intakes
	minimumDryMass += 0.01 //Hinge
	//minimumDryMass += 0.001 //Occlude
	minimumDryMass += 0.08 //Control
	minimumDryMass += 0.24 * fcArrays //Fuel Cell test
	minimumDryMass += 0.025 //Heat shielding

	mk0s := 0.0
	if lf > 0.5 {
		mk0s = math.Ceil(4.0 * (lf - 0.5))
		minimumDryMass += mk0s * 0.025
	}

	minXe := minXeCapacity(totalXe)
	minimumDryMass += minXe / 3.0 //Xenon tanks

	minLFOX := minLFOXCapacity(totalLFOX)
	fmt.Printf("Xenon Capacity: %ft (%f units)\nLF+OX Mixture Capacity: %ft (%f units)\nLF Capacity: %ft (%f units)\n\n", minXe, 10000.0 * minXe, minLFOX, 200.0 * minLFOX, 0.5 + 0.25 * mk0s, 100.0 + 50.0 * mk0s)
	//fmt.Print(minXe, minLFOX, 0.5 + 0.25 * mk0s, "\n")

	minimumDryMass += minLFOX / 8.0 //LFOX tanks

	if dryMass >= minimumDryMass {
		fmt.Printf("Calculated Dry Mass: %ft (excess of %dkg / %fm/s)\n", minimumDryMass, int(math.Round(1000.0*(dryMass - minimumDryMass))), 4200.0*g*math.Log(dryMass/minimumDryMass))
		fmt.Printf("Predicted Minimum Mass: %ft\n", startingMass * minimumDryMass/dryMass)
	} else {
		fmt.Printf("Calculated Dry Mass: %ft (deficit of %dkg)\n", minimumDryMass, int(math.Round(1000.0*(minimumDryMass - dryMass))))
	}
}