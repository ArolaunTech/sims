import matplotlib.pyplot as plt

startcoolantTemp = 44
startthermiteTemp = 2700
transferFraction = 0.5 #Between 0.5 and 1

res = 24
coolantTemps = [0 for i in range(res)]
thermiteTemps = [0 for i in range(res)]

for i in range(4000):
	#Converge
	newCoolantTemps = [coolantTemps[0]]
	for i in range(1,res):
		newCoolantTemps.append(transferFraction * coolantTemps[i-1] + (1-transferFraction) * thermiteTemps[i])
	newThermiteTemps = [0 for i in range(res)]
	newThermiteTemps[-1] = thermiteTemps[-1]
	for i in range(res-1):
		newThermiteTemps[i] = transferFraction*thermiteTemps[i+1] + (1-transferFraction)*coolantTemps[i]

	coolantTemps = newCoolantTemps
	thermiteTemps = newThermiteTemps
	#Shift
	thermiteTemps[:-1] = thermiteTemps[1:]
	coolantTemps[1:] = coolantTemps[:-1]
	#Push
	coolantTemps[0] = startcoolantTemp
	thermiteTemps[-1] = startthermiteTemp

print(coolantTemps[0], thermiteTemps[0])
print(coolantTemps[-1], thermiteTemps[-1])
print('\n')
print(thermiteTemps[0] - coolantTemps[0])

#plt.plot(coolantTemps)
#plt.plot(thermiteTemps)
#plt.show()