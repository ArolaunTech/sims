import numpy as np

generator = np.random.default_rng()

g0 = 9.81

R = 200000
initalt = 1
gsurf = 1.63

mu = R * R * gsurf

thrust = 20000
initmass = 10000
isp = 320

dt = 0.1

def sign(x):
	if x > 0:
		return 1
	if x < 0:
		return -1
	return 0

def simulate(ipx, ipy, idpx, idpy, log=False):
	x = 0
	y = R + initalt
	vx = 0
	vy = 0

	px = ipx
	py = ipy
	dpx = idpx
	dpy = idpy

	mass = initmass
	time = 0

	for i in range(10000):
		dist = np.sqrt(x * x + y * y)
		grav = mu / dist / dist
		gx = -grav * x / dist
		gy = -grav * y / dist

		dgxdx = mu * (2 * x * x - y * y) / dist ** 5
		dgxdy = 3 * mu * x * y / dist ** 5
		dgydx = dgxdy
		dgydy = mu * (2 * y * y - x * x) / dist ** 5

		ddpx = dgxdx * px + dgxdy * py
		ddpy = dgydx * px + dgydy * py

		dpx += ddpx * dt
		dpy += ddpy * dt
		px += dpx * dt
		py += dpy * dt

		pmag = np.sqrt(px * px + py * py)
		if pmag > 1:
			accel = thrust / mass

			gx += accel * px / pmag
			gy += accel * py / pmag
			mass -= thrust / isp / g0 * dt

			if mass < 0:
				return False, 0

		vx += gx * dt
		vy += gy * dt
		x += vx * dt
		y += vy * dt

		h = y * vx - x * vy
		a = 1 / (2 / dist - (vx * vx + vy * vy) / mu)
		e = np.sqrt(1 - h * h / mu / a)
		vh = h / dist
		vv = (x * vx + y * vy) / dist

		if a < 0:
			return False, 0
		if a * (1 - e) > R:
			break
		if dist < R:
			return False, 0

		dv = isp * g0 * np.log(initmass / mass)

		pv = (x * px + y * py) / dist
		ph = (y * px - x * py) / dist

		if log and i % 100 == 0:
			print(f"Time {round(time, 1)}s:\nAlt: {round(dist) - R}m\tvh: {round(vh, 2)} m/s\tvv: {round(vv, 2)} m/s\tAngle: {round(np.atan2(pv, ph) * 180 / np.pi, 3)} deg")
		time += dt
	if log:
		print(f"Time {round(time, 1)}s:\nAlt: {round(dist) - R}m\tvh: {round(vh, 2)} m/s\tvv: {round(vv, 2)} m/s\tAngle: {round(np.atan2(pv, ph) * 180 / np.pi, 3)} deg")

	if a * (1 - e) <= R:
		return False, 0

	realdv = dv + np.sqrt(mu / R) - np.sqrt(mu * (2 / R - 1 / a))

	if log:
		print(f"dv used: {round(dv, 3)} m/s\tapoapsis: {round(a * (1 + e) - R)} m\tfinal mass: {round(mass, 3)}/{initmass} kg")

	#print(round(vh, 3), round(vv, 3), round(dist, 3), round(initmass * np.exp(-realdv / isp / g0), 3), round(realdv, 3))
	return True, dv

def constalt():
	vel = 0

	mass = initmass

	grav = mu / (R + initalt) ** 2

	while vel < np.sqrt(mu / (R + initalt)):
		vaccel = grav - vel * vel / (R + initalt)
		haccel = np.sqrt((thrust / mass) ** 2 - vaccel * vaccel)

		mass -= thrust / isp / g0 * dt
		vel += haccel * dt

	return isp * g0 * np.log(initmass / mass), mass

def gravturn(initangle, log=False):
	x = 0
	y = R + initalt
	vy = 1e-3
	vx = vy * initangle

	mass = initmass

	for i in range(10000):
		dist = np.sqrt(x * x + y * y)

		if dist < R:
			return -1, vv

		vel = np.sqrt(vx * vx + vy * vy)

		grav = mu / dist / dist
		gx = -x * grav / dist
		gy = -y * grav / dist

		accel = thrust / mass
		ax = accel * vx / vel
		ay = accel * vy / vel

		vx += (gx + ax) * dt
		vy += (gy + ay) * dt
		x += vx * dt
		y += vy * dt
		mass -= thrust / isp / g0 * dt

		h = y * vx - x * vy
		a = 1 / (2 / dist - (vx * vx + vy * vy) / mu)
		e = np.sqrt(1 - h * h / mu / a)
		vh = h / dist
		vv = (x * vx + y * vy) / dist

		if a < 0 and vv > 0:
			return -2, vv
		if a * (1 - e) > R:
			break

		if log and i % 100 == 0:
			print(dist, vh, vv, mass)

	return mass, vv

print(constalt())

lo = -21
hi = 0
for i in range(0):
	#print(lo, hi)

	mid = (lo + hi) / 2
	res, vv = gravturn(np.exp(mid))

	#print(res, vv, mid)

	if res == -1:
		hi = mid
		continue
	if res == -2:
		lo = mid
		continue
	if vv < -1e-6:
		hi = mid
		continue
	if vv > 1e-6:
		lo = mid
		continue

	break

#print(np.exp(mid))
#print(gravturn(np.exp(mid), True))

bresult = False
bscore = 0

while not bresult:
	bpx = generator.normal(scale=2)
	bpy = generator.normal(scale=2)
	bdpx = generator.normal(scale=0.02)
	bdpy = generator.normal(scale=0.02)

	bresult, bscore = simulate(bpx, bpy, bdpx, bdpy)

for i in range(100000):
	npx = bpx + generator.normal(scale=0.0002)
	npy = bpy + generator.normal(scale=0.0002)
	ndpx = bdpx + generator.normal(scale=0.000002)
	ndpy = bdpy + generator.normal(scale=0.000002)

	nresult, nscore = simulate(npx, npy, ndpx, ndpy)

	if not nresult:
		continue
	if nscore >= bscore:
		continue

	bpx = npx
	bpy = npy
	bdpx = ndpx
	bdpy = ndpy

	bscore = nscore

	print(bpx, bpy, bdpx, bdpy, bscore, initmass * np.exp(-bscore / isp / g0))

	simulate(bpx, bpy, bdpx, bdpy, True)
	print(i)