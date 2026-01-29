import rapiermaxspeed
import math

g0 = 9.81
pi = 3.14159265358979323846264338

R = 600000
stdgp = g0 * R * R
rotperiod = 21549.425
rotvel = 2 * pi * R / rotperiod

coeffs = [0.41]

mass = 2.85
thrust = 40
isp = 320
fc = thrust / isp / g0

t = 0
dt = 1

x = 0
y = R

vx = 710 * math.cos(coeffs[0]) + rotvel
vy = 710 * math.sin(coeffs[0])

for i in range(10000):
	energy = (vx * vx + vy * vy) / 2 - stdgp / math.sqrt(x * x + y * y)
	sma = -stdgp / 2 / energy
	h = y * vx - x * vy
	eccen = math.sqrt(1 - (h * h) / stdgp / sma)

	if sma * (1 + eccen) > R + 145000:
		break

	print(round(t, 1), round(math.sqrt(x * x + y * y) - R), round(mass, 3), round(sma * (1 + eccen) - R), round(sma * (1 - eccen) - R), round(eccen, 3))

	sqdist = x * x + y * y
	dist = math.sqrt(sqdist)

	grav = stdgp / sqdist
	gx = -x * grav / dist
	gy = -y * grav / dist

	# Calculate thrust offset
	sx = 2 * pi * y / rotperiod
	sy = -2 * pi * x / rotperiod

	vh = (y * vx - x * vy) / dist - 2 * pi * dist / rotperiod
	vv = (x * vx + y * vy) / dist - grav * dt

	a = sum([coeffs[j] * t ** j for j in range(len(coeffs))])
	ca = math.cos(a)
	sa = math.sin(a)

	sspeed = math.sqrt(vh * vh + vv * vv)
	temp = rapiermaxspeed.temperatureAlt(dist - R)
	sound = math.sqrt(rapiermaxspeed.KERBIN_ADIABATIC * rapiermaxspeed.R * temp / rapiermaxspeed.KERBIN_MOLARMASS)

	ld = rapiermaxspeed.ldratio.evaluate(sspeed / sound)

	p0x = (-ld * vv - vh) / (sspeed * math.sqrt(ld * ld + 1))
	p0y = (ld * vh - vv) / (sspeed * math.sqrt(ld * ld + 1))
	p0x2 = (ld * vv - vh) / (sspeed * math.sqrt(ld * ld + 1))
	p0y2 = (-ld * vh - vv) / (sspeed * math.sqrt(ld * ld + 1))

	T = thrust / mass * dt

	c0 = p0x * ca + p0y * sa
	c1 = -p0x * sa + p0y * ca
	c2 = -vh * sa + vv * ca
	c3 = vh * ca + vv * sa
	c02 = p0x2 * ca + p0y2 * sa
	c12 = -p0x2 * sa + p0y2 * ca

	L = (c0 * T - c2) / c1
	L2 = (c02 * T - c2) / c12
	O = c3 + c0 * L + math.sqrt(T * T - (c2 + c1 * L) ** 2)
	O2 = c3 + c02 * L2 + math.sqrt(T * T - (c2 + c12 * L2) ** 2)

	Oreal = 0
	if L > 0:
		Oreal = O
	if L2 > 0 and O2 > Oreal:
		Oreal = O2

	newvx = (Oreal * math.cos(a) + 2 * pi * dist / rotperiod) * y / dist + Oreal * math.sin(a) * x / dist
	newvy = -(Oreal * math.cos(a) + 2 * pi * dist / rotperiod) * x / dist + Oreal * math.sin(a) * y / dist

	t += dt
	vx = newvx
	vy = newvy
	mass -= fc * dt
	x += vx * dt
	y += vy * dt