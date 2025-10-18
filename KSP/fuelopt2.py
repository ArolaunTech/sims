import math

tankinfo = {
	"None": (0, 0, 0),
	"Dumpling": (0.01375, 0.11, 0.12375),
	"Oscar": (0.025, 0.2, 0.225),
	"Baguette": (0.03375, 0.27, 0.30375),
	"Donut": (0.0375, 0.3, 0.3375),
	"FL-T100": (0.0625, 0.5, 0.5625)
}

g0 = 9.80665

payload = 0.05 + 0.07 + 0.02 + 0.02 + 0.005
dv = 500 + 900 + 1400 + 500
decoupler = 0.01

isp = 301

def calcdv(tanks):
	out = 0

	dry = payload
	wet = payload

	for tank in tanks:
		dry += tankinfo[tank][0]
		wet += tankinfo[tank][2]

		out += isp * g0 * math.log(wet / dry)

		wet += decoupler
		dry = wet

	return out

def calcscore(tanks):
	return payload + sum([tankinfo[tank][2] for tank in tanks]) + decoupler * (len(tanks) - 1)

for i in range(216):
	tanks = [list(tankinfo.keys())[i % 6], list(tankinfo.keys())[(i // 6) % 6], list(tankinfo.keys())[(i // 36) % 6]]
	tankdv = calcdv(tanks)

	if tankdv < dv:
		continue

	print(tanks, tankdv, calcscore(tanks))