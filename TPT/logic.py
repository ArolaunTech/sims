def f(a):
	return [a[i] for i in range(len(a)) if (i==0) or a[i] != a[i-1]]

def addNew(value, val, i):
	if informationCondition(value):
		if tuple(value) not in already_found:
			already_found.add(tuple(value))
			newVals.append(value)
			if finalCondition(value):
				print(value, val, i)
			elif (1 in value) or (3 in value):
				print('Close!', value, val, i)

def finalCondition(value):
	return (1 in value) and (3 in value)

def informationCondition(value):
	return (1 in value) or (2 in value)

oldVals = [[1,3]]
print(oldVals)
already_found = {(1,3)}
numSearched = 0
maxBitValue = 3

for i in range(20):
	newVals = []
	for val in oldVals:
		for j in range(len(val)):
			addNew(f(val[:j] + val[j+1:]), val, i)
			#if val[j] < 3:
			#	addNew(f(val[:j] + [maxBitValue - val[j]] + val[j:]), val, i)
			if (j > 0) and (val[j-1] & (maxBitValue - val[j]) > 0):
				addNew(f(val[:j] + [val[j-1] & (maxBitValue - val[j])] + val[j:]), val, i)
				#addNew(val[:j] + [val[j-1] & (maxBitValue - val[j])] + val[j:])
	oldVals = newVals
	numSearched += len(newVals)
	print(numSearched)
	#print(newVals)
	#print(oldVals, newVals)