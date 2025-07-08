import tensorflow as tf
import keras
import numpy as np
from copy import deepcopy

from utils import *

#Constants
timestep = 0.1

#Engine sorting
engines_sorted = sorted([engine for engine in engine_info], key=lambda x: engine_info[x][0] * (1 + 1 / fuel_tank_info[engine_info[x][1]][1]), reverse=True)

#Simulation
def simulate_ascent(grav_parameter, radius, rotational_period, start_alt, fuel_reverse, craft_mass, craft_engines, agent, predictor):
	x = 0
	y = radius + start_alt
	vx = 2 * np.pi * y / rotational_period
	vy = 0
	mass = craft_mass

	total_thrust = sum([craft_engines[engine] * engine_info[engine][2] for engine in craft_engines])

	model_inputs = [craft_engines[engine] for engine in engine_info]
	model_inputs.append(grav_parameter / radius / radius)
	model_inputs.append(np.sqrt(grav_parameter / radius))
	model_inputs.append(radius)
	model_inputs.append(radius + start_alt)
	model_inputs.append(vx)
	model_inputs.append(vy)
	model_inputs.append(mass)

	engines = [engine for engine in engine_info]

	agent_inputs = []
	agent_outputs = []

	periapsis = 0
	time = 0
	while periapsis < radius:
		if time > 10000:
			break

		#Calculate acceleration
		dist = np.sqrt(x * x + y * y)
		if dist < radius:
			break

		vh = (-y * vx + x * vy) / dist
		vv = (x * vx + y * vy) / dist

		model_inputs[-4] = dist
		model_inputs[-3] = vh
		model_inputs[-2] = vv
		model_inputs[-1] = mass

		agent_inputs.append(deepcopy(model_inputs))

		throttle, th, tv = tuple(agent.predict(np.array(model_inputs).reshape((-1, 33)), verbose=0)[0])
		agent_outputs.append([throttle, th, tv])

		if throttle > 10:
			thrust = total_thrust
		elif throttle < -10:
			thrust = 0
		else:
			thrust = total_thrust / (1 + np.exp(-throttle))

		accel = thrust / mass
		ah = accel * th / np.sqrt(th * th + tv * tv)
		av = accel * tv / np.sqrt(th * th + tv * tv)

		ax = (-y * ah + x * av) / dist
		ay = (x * ah + y * av) / dist

		tick_consumption = {fuel: 0 for fuel in fuel_tank_info}
		for engine in engines_sorted:
			if engine not in craft_engines:
				continue

			group_thrust = craft_engines[engine] * engine_info[engine][2]

			if group_thrust > thrust:
				tick_consumption[engine_info[engine][1]] += thrust / engine_info[engine][0] / 9.80665
				break

			thrust -= group_thrust
			tick_consumption[engine_info[engine][1]] += group_thrust / engine_info[engine][0] / 9.80665

		#Gravity
		grav = grav_parameter / dist / dist
		gx = -x * grav / dist
		gy = -y * grav / dist

		consumption_mult = -1 if fuel_reverse else 1
		for consumption in tick_consumption.values():
			mass -= consumption_mult * consumption * timestep

		vx += (ax + gx) * timestep
		vy += (ay + gy) * timestep

		x += vx * timestep
		y += vy * timestep

		time += timestep

		sma = 1 / (2 / dist - (vx * vx + vy * vy))
		eccentricity = np.sqrt(1 - (dist * dist * vh * vh)/ sma / grav_parameter)

		periapsis = sma * (1 - eccentricity)

	if periapsis < radius:
		return [periapsis/radius, agent_inputs, agent_outputs]
	return [1, agent_inputs, agent_outputs]

def simulate_landing(grav_parameter, radius, rotational_period, land_alt, orbit_alt, craft_mass, craft_engines, agent, predictor):
	low = 0
	high = craft_mass

	while high - low > 1e-3:
		middle = (low + high) / 2

		result = simulate_ascent(grav_parameter, radius, rotational_period, land_alt, True, middle, craft_engines, agent, predictor)

		print([a + b for a, b in zip(result[1], result[2])])

		low = middle


#Model definitions
def get_agent():
	inputs = keras.Input(shape=(len(engine_info) + 7,), name="state")
	x1 = keras.layers.Dense(32, activation="relu")(inputs)
	x2 = keras.layers.Dense(32, activation="relu")(x1)
	outputs = keras.layers.Dense(3, name="action")(x2)

	model = keras.Model(inputs=inputs, outputs=outputs)
	return model
	
def get_predictor():
	inputs = keras.Input(shape=(len(engine_info) + 10,), name="stateaction")
	x1 = keras.layers.Dense(32, activation="relu")(inputs)
	x2 = keras.layers.Dense(32, activation="relu")(x1)
	outputs = keras.layers.Dense(1, name="reward")(x2)

	model = keras.Model(inputs=inputs, outputs=outputs)
	return model

#Get networks
agent = get_agent()
predictor = get_predictor()

agent.summary()

#Training classes
optimizer = keras.optimizers.SGD(learning_rate=1e-3)
loss = keras.losses.MeanSquaredError()

#Train
for i in range(1000):
	radius = 2e5
	grav_parameter = 65138397520.7807
	rotational_period = 138984.376574476

	craft_engines = {engine: 0 for engine in engine_info}
	craft_engines["Vector"] = 1

	simulate_landing(grav_parameter, radius, rotational_period, 1000, 10000, 1e5, craft_engines, agent, predictor)