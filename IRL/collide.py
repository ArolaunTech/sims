sim = [
	0, 0, 0, 
	0, 1, 0,
	0, 0, 0
]
si = sum(sim)

df = 0.001
diffuse = [
	[5,7,4,8],
	[3,7,4,6],
	[1,3,4,0],
	[1,5,4,2],
	[0,4,1,3],
	[2,4,1,5],
	[6,4,3,7],
	[8,4,5,7]
]
T=1

mdiff = 1e-10

for i in range(9):
	sim[i] += T

for i in range(100):
	new = [0, 0, 0, 0, 0, 0, 0, 0, 0]
	for i in range(9):
		new[i] = sim[i]
	#Do stuff
	for reac in diffuse:
		am = sim[reac[0]] * sim[reac[1]]
		am *= df
		new[reac[2]] += am
		new[reac[3]] += am
		new[reac[0]] -= am
		new[reac[1]] -= am
	diff = 0
	for i in range(9):
		diff += abs(sim[i] - new[i])
	if diff < mdiff:
		break
	sim = [0, 0, 0, 0, 0, 0, 0, 0, 0]
	for i in range(9):
		sim[i] = new[i]
print(sim)
s = sum(sim)
for i in range(4):
	mi = min(sim[i], sim[8 - i])
	sim[i] -= mi
	sim[8 - i] -= mi
s = sum(sim)
sim[4] = max(0, sim[4] - s + si)
s = sum(sim)
for i in range(9):
	sim[i] *= si / s
for i in range(9):
	sim[i] = round(sim[i], 3)
print(sim)
print(sim[2] + sim[5] + sim[8] - sim[0] - sim[3] - sim[6])