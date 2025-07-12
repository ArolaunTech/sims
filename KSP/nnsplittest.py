import os
os.environ["KERAS_BACKEND"] = "tensorflow"

import tensorflow as tf
import keras
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

g0 = 9.80665
timestep = 0.1

generator = np.random.default_rng()

planetRadius = 300000
planetGrav = planetRadius ** 2 * 0.275 * g0

def simulate(initDist, initHv, initVv, ax, ay, isp, time):
	x = initDist
	y = 0
	vx = initVv
	vy = initHv
	currTime = 0

	accel = np.sqrt(ax * ax + ay * ay)
	accelnormalized = accel / isp / g0
	dvnormalized = 0

	states = []

	while currTime < time:
		distSquared = x * x + y * y
		dist = np.sqrt(distSquared)

		vh = (-y * vx + x * vy) / dist
		vv = (x * vx + y * vy) / dist

		grav = planetGrav / distSquared
		gx = -x * grav / dist
		gy = -y * grav / dist

		vx += (ax + gx) * timestep
		vy += (ay + gy) * timestep
		x += vx * timestep
		y += vy * timestep
		currTime += timestep
		dvnormalized += accelnormalized * timestep

		states.append([dist * 1e-6, vh * 1e-3, vv * 1e-3, accel * np.exp(dvnormalized)])

	return states

def split(dist1, dist2, hv1, hv2, vv1, vv2, accel, isp, model, level):
	if level <= 0:
		modelinput = np.array([
			dist1 * 1e-6, 
			hv1 * 1e-3, 
			vv1 * 1e-3, 
			np.log(accel), 
			dist2 * 1e-6,
			hv2 * 1e-3,
			vv2 * 1e-3
		]).reshape((-1,7))

		modeloutput = model(modelinput)

		#print(modeloutput)

		return np.array(modeloutput)

	dsplits, hsplits, vsplits = [], [], []
	for i in range(level + 2):
		frac = i / (level + 1)
		dsplits.append(frac * dist2 + (1 - frac) * dist1)
		hsplits.append(frac * hv2   + (1 - frac) * hv1)
		vsplits.append(frac * vv2   + (1 - frac) * vv1)

	newaccel = accel
	for i in range(level + 1):
		prevaccel = newaccel
		newaccel = split(dsplits[i], dsplits[i + 1], hsplits[i], hsplits[i + 1], vsplits[i], vsplits[i + 1], newaccel, isp, model, 0)[0][0]

		if newaccel < prevaccel:
			newaccel = prevaccel

	for i in range(100):
		newdsplits, newhsplits, newvsplits = deepcopy((dsplits, hsplits, vsplits))

		for j in range(1, level + 1):
			newdsplits[j] += generator.normal(loc=0, scale=1000)
			newhsplits[j] += generator.normal(loc=0, scale=100)
			newvsplits[j] += generator.normal(loc=0, scale=50)

		newnewaccel = accel
		for i in range(level + 1):
			prevaccel = newnewaccel
			newnewaccel = split(newdsplits[i], newdsplits[i + 1], newhsplits[i], newhsplits[i + 1], newvsplits[i], newvsplits[i + 1], newnewaccel, isp, model, 0)[0][0]

			if newnewaccel < prevaccel:
				newnewaccel = prevaccel

		if newnewaccel < newaccel:
			newaccel = newnewaccel

			dsplits, hsplits, vsplits = deepcopy((newdsplits, newhsplits, newvsplits))

	return np.array([[newaccel]])

#NN
def nnloss(actual, predicted):
	diff = predicted - actual
	diffnormalized = tf.math.maximum(tf.math.minimum(diff, 40), -40)
	diffclamped = tf.math.maximum(tf.math.minimum(diff, 5), -5)

	loss = tf.math.exp(-1 * diffnormalized * diffnormalized) * diffclamped + 4e-5 * diff * diff + 0.5

	return tf.math.reduce_mean(25 * loss)

def getModel():
	inputs = keras.Input(shape=(7,), name="boundaries")
	x1 = keras.layers.Dense(16, activation="relu")(inputs)
	x2 = keras.layers.Dense(16, activation="relu")(x1)
	outputs1 = keras.layers.Dense(1, name="prediction")(x2)
	#outputs2 = keras.layers.Dense(1, activation="relu", name="std-dev")(x2)

	model = keras.Model(inputs=inputs, outputs=outputs1)
	return model

model = getModel()
optimizer = keras.optimizers.Adam(learning_rate=1e-3)
lossFn = nnloss

model.summary()

#Train
losses = []
for i in range(1000):
	if generator.random() < 0.3:
		#Split
		initDist = planetRadius + generator.uniform(low=0, high=10000)
		initHv = generator.normal(loc=0, scale=500)
		initVv = generator.normal(loc=0, scale=20)

		finalDist = planetRadius + generator.uniform(low=0, high=10000)
		finalHv = generator.normal(loc=0, scale=500)
		finalVv = generator.normal(loc=0, scale=20)

		accel = abs(generator.normal(loc=0, scale=5))
		if accel < 1:
			accel = 1

		level = generator.integers(16)

		inputs = np.array([
			initDist,
			initHv,
			initVv,
			np.log(accel),
			finalDist,
			finalHv,
			finalVv
		])

		outputs = split(initDist, finalDist, initHv, finalHv, initVv, finalVv, accel, 320, model, level)

		with tf.GradientTape() as tape:
			modeloutputs = model(inputs.reshape((-1,7)), training=True)

			loss = nnloss(outputs, modeloutputs)
		grads = tape.gradient(loss, model.trainable_weights)
		optimizer.apply(grads, model.trainable_weights)

		losses.append(np.array(loss))
		print(i, 0, np.array(loss))
		continue

	#Simulation
	initDist = planetRadius + generator.uniform(low=0, high=10000)
	initHv = generator.normal(loc=0, scale=500)
	initVv = generator.normal(loc=0, scale=20)

	ax = generator.normal(loc=0, scale=5)
	ay = generator.normal(loc=0, scale=5)
	while ax * ax + ay * ay < 1:
		ax = generator.normal(loc=0, scale=5)
		ay = generator.normal(loc=0, scale=5)

	timelength = np.exp(generator.normal(loc=5, scale=1))

	datapoints = simulate(initDist, initHv, initVv, ax, ay, 320, timelength)
	numdatapoints = len(datapoints)

	for j in range(10):
		inputs = []
		outputs = []
		for k in range(32):
			start = generator.integers(numdatapoints)
			end = generator.integers(numdatapoints)

			while start == end:
				start = generator.integers(numdatapoints)
				end = generator.integers(numdatapoints)

			if start > end:
				temp = start
				start = end
				end = temp

			startpoint = datapoints[start]
			endpoint = datapoints[end]

			inputs.append(np.array(startpoint[:3] + [np.log(startpoint[3])] + endpoint[:3]))
			outputs.append(np.log(endpoint[3]))

		inputs = np.array(inputs)
		outputs = np.array(outputs).reshape(-1,1)

		with tf.GradientTape() as tape:
			modeloutputs = model(inputs.reshape((-1,7)), training=True)

			loss = nnloss(outputs, modeloutputs)
		grads = tape.gradient(loss, model.trainable_weights)
		optimizer.apply(grads, model.trainable_weights)

		#print(i, j, np.array(loss), inputs, outputs, modeloutputs)

	losses.append(np.array(loss))
	print(i, 1, np.array(loss))

plt.plot(losses)
plt.yscale('log')
plt.show()