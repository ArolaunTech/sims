import os
os.environ["KERAS_BACKEND"] = "tensorflow"

import tensorflow as tf
import keras
import tensorflow.keras.backend as backend
import numpy as np
import matplotlib.pyplot as plt

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

#NN
def nnloss(actual, predicted):
	mean = predicted[0]
	stddev = predicted[1]

	stddev = backend.maximum(stddev, 0.1)

	negativelogprobabilitydensity = 0.5 * (np.log(2 * np.pi) + 2 * stddev) + backend.square(actual - mean) / (2 * stddev * stddev)

	return backend.mean(negativelogprobabilitydensity)

def getModel():
	inputs = keras.Input(shape=(7,), name="boundaries")
	x1 = keras.layers.Dense(16, activation="relu")(inputs)
	x2 = keras.layers.Dense(16, activation="relu")(x1)
	outputs1 = keras.layers.Dense(1, name="mean")(x2)
	outputs2 = keras.layers.Dense(1, activation="relu", name="std-dev")(x2)

	model = keras.Model(inputs=inputs, outputs=[outputs1, outputs2])
	return model

model = getModel()
optimizer = keras.optimizers.Adam(learning_rate=1e-3)
lossFn = nnloss

model.summary()

#Train
losses = []
for i in range(100):
	initDist = planetRadius + generator.uniform(low=0, high=10000)
	initHv = generator.normal(loc=0, scale=500)
	initVv = generator.normal(loc=0, scale=20)

	ax = generator.normal(loc=0, scale=5)
	ay = generator.normal(loc=0, scale=5)

	timelength = np.exp(generator.normal(loc=5, scale=1))

	datapoints = simulate(initDist, initHv, initVv, ax, ay, 320, timelength)
	numdatapoints = len(datapoints)

	for j in range(100):
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

			inputs.append(np.array(startpoint + endpoint[:3]))
			outputs.append(endpoint[3])

		inputs = np.array(inputs)
		outputs = np.array(outputs).reshape(-1,1)

		with tf.GradientTape() as tape:
			modeloutputs = model(inputs.reshape((-1,7)), training=True)

			loss = nnloss(outputs, modeloutputs)
		grads = tape.gradient(loss, model.trainable_weights)
		optimizer.apply(grads, model.trainable_weights)

		#print(i, j, np.array(loss), inputs, outputs, modeloutputs)

	losses.append(np.array(loss))
	print(i, np.array(loss))

plt.plot(losses)
plt.yscale('log')
plt.show()