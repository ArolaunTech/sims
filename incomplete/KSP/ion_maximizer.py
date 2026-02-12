import matplotlib.pyplot as plt
import numpy as np

gr = (np.sqrt(5) + 1)/2
ionEC = 600/(7 * 9.81)

def calcStarvationEngine(throttle, ec):
	#Returns thrust (kN), EC (zaps/s), Xenon (tons/s)
	if ec <= 0: #No EC
		return 0, 0, 0
	if throttle <= 0:
		return 0, 0, 0
	if throttle > 1:
		return calcStarvationEngine(1, ec)

	#Here we know that EC > 0 and 0 < throttle < 1
	demandXenon = 2 * throttle/(4200 * 9.81)
	demandEC = 180000 * demandXenon

	if ec > demandEC:
		return 2 * throttle, demandEC, demandXenon

	#Ion starvation logic
	satisfaction = ec/demandEC
	runFrac = 1
	if satisfaction < 0.1:
		runFrac = 10 * satisfaction
		satisfaction = 0.1

	demandXenon *= runFrac
	ispMult = 0.5 + 0.5 * satisfaction

	return 2 * throttle * ispMult * runFrac, ec, demandXenon

def cninter(ec, nEngines, tfirst):
	tf, ef, xf = calcStarvationEngine(tfirst, ec)
	ts, es, xs, pop = calcNEngines(ec - ef, nEngines - 1)

	return tf + ts, ef + es, xf + xs, pop

def calcNEngines(ec, nEngines):
	if ec <= 0:
		return 0, 0, 0, 1
	if ec >= nEngines * 600/(7 * 9.81):
		return 2 * nEngines, nEngines * 600/(7 * 9.81), nEngines * 2/(4200 * 9.81), 1
	if nEngines == 1:
		return calcStarvationEngine(1, ec) + tuple([1])

	a, b, c, d = 0, 2 - gr, gr - 1, 1
	tar, ear, xar, par = cninter(ec, nEngines, a)
	if par < 1e-4:
		return tar, ear, xar, 0
	tbr, ebr, xbr, pbr = cninter(ec, nEngines, b)
	tcr, ecr, xcr, pcr = cninter(ec, nEngines, c)
	tdr, edr, xdr, pdr = cninter(ec, nEngines, d)

	for i in range(10):
		#print(tar, tbr, tcr, tdr)
		#print(a, b, c, d)
		if max(tar, tbr) > max(tcr, tdr):
			tdr, edr, xdr, d = tcr, ecr, xcr, c
			tcr, ecr, xcr, c = tbr, ebr, xbr, b
			b = d * (2 - gr) + a * (gr - 1)
			tbr, ebr, xbr, pbr = cninter(ec, nEngines, b)
		else:
			tar, ear, xar, a = tbr, ebr, xbr, b
			tbr, ebr, xbr, b = tcr, ecr, xcr, c
			c = d * (gr - 1) + a * (2 - gr)
			tcr, ecr, xcr, pcr = cninter(ec, nEngines, c)

	
	#print(nEngines)
	return tar, ear, xar, a

def calcNEnginesOpt(ec, nEngines):
	if ec <= 0:
		return 0, 0, 0
	if ec >= nEngines * ionEC:
		return 2 * nEngines, nEngines * ionEC, nEngines * 2/(4200 * 9.81)

	maxSU = min(ec, ionEC * (nEngines - 1))
	minSU = max(0, ec - ionEC)

	lU = ec - maxSU
	uU = ec - minSU

	zU = 0.1 * ionEC
	if uU < 0.1 * ionEC:
		zU = uU
	if lU > 0.1 * ionEC:
		zU = lU

	res = calcStarvationEngine(1, zU)
	bestT = 2 * (ec - zU)/ionEC + res[0]
	bestXe = 2 * (ec - zU)/(ionEC * 4200 * 9.81) + res[2]
	#print(ec, bestT, bestXe)

	#bestT, bestEC, bestXe = 2 * wholeSatisfied, ec, 2 * wholeSatisfied/(4200*9.81)
	return bestT, ec, bestXe

xs = np.linspace(0.001, 70, num=100000)
txs = []
exs = []
xxs = []
for x in xs:
	tx, ex, xx = calcNEnginesOpt(x, 7)
	print(tx, ex, xx)
	txs.append(tx)
	exs.append(tx/(9.81*xx))
	xxs.append(xx)

plt.plot(xs, exs)
plt.xlabel("EC consumption (ec/s)")
plt.ylabel("Specific Impulse (seconds)")
plt.title("Isp of 7 EC-starved Ion Engines Optimized for Maximum Thrust")
plt.show()

def calcTwoEngines(throttle1, ec):
	t1, e1, x1 = calcStarvationEngine(throttle1, ec)
	t2, e2, x2 = calcStarvationEngine(1, ec - e1)
	return t1 + t2, e1 + e2, x1 + x2

ec = np.linspace(0, 20, num=100)
throttle = np.linspace(0, 1, num=100)

ec, throttle = np.meshgrid(ec, throttle)

starveCalcVec = np.vectorize(calcStarvationEngine, otypes=[np.float64, np.float64, np.float64])
calcTwoVec = np.vectorize(calcTwoEngines, otypes=[np.float64, np.float64, np.float64])

thrust, ecConsump, xenonConsump = calcTwoVec(throttle, ec)
maxThrust = np.max(thrust, axis=0)
print(list(maxThrust))

plt.imshow(thrust, cmap='inferno')
#plt.show()