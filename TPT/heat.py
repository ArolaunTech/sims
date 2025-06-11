import random

elements = {
	"NONE": (0, 0, False),
	"CRMC": (5164 + 273.15, 35, True),
	"TUNG": (3695.0, 251, True),
	"PTNM": (1768 + 273.15, 251, True),
	"GLOW": (1e+6, 44, True),
	"PLSM": (1e+6, 5, True)
}

ROWS = 50
COLS = 5

elementMap = [["GLOW", "CRMC", "PLSM", "PLSM", "PLSM", "PLSM"] for i in range(ROWS)]

ROWS = len(elementMap)
COLS = len(elementMap[0])

maxtemp = [[elements[elementMap[i][j]][0] for j in range(COLS)] for i in range(ROWS)]
conductivity = [[elements[elementMap[i][j]][1] for j in range(COLS)] for i in range(ROWS)]
existence = [[elements[elementMap[i][j]][2] for j in range(COLS)] for i in range(ROWS)]

failures = 0

for sim in range(10):
	#Init
	temp = [[9999 if elementMap[i][j] == "PLSM" else 22 + 273.15 for j in range(COLS)] for i in range(ROWS)]
	failure = False

	for iteration in range(100):
		#Iterate over everything but plasma
		for i in range(ROWS):
			for rj in range(COLS):
				j = (rj + 1)%COLS
				if elementMap[i][j] == "PLSM" or elementMap[i][j] == "NONE":
					continue
				if random.random() > conductivity[i][j]/250:
					continue
				hcount = 1
				cheat = temp[i][j]
				surround = [(i, j)]
				for di in range(-1, 2):
					for dj in range(-1, 2):
						ni = i + di
						nj = j + dj
						if ni < 0 or ni >= ROWS:
							continue
						if nj < 0 or nj >= COLS:
							continue
						if elementMap[ni][nj] == "NONE":
							continue
						hcount += 1
						cheat += temp[ni][nj]
						surround.append((ni, nj))
				cheat /= hcount
				for x, y in surround:
					temp[x][y] = cheat

		for i in range(ROWS):
			for j in range(COLS):
				if elementMap[i][j] != "PLSM":
					continue
				if random.random() > conductivity[i][j]/250:
					continue
				hcount = 1
				cheat = temp[i][j]
				surround = [(i, j)]
				for di in range(-1, 2):
					for dj in range(-1, 2):
						ni = i + di
						nj = j + dj
						if ni < 0 or ni >= ROWS:
							continue
						if nj < 0 or nj >= COLS:
							continue
						if elementMap[ni][nj] == "NONE":
							continue
						hcount += 1
						cheat += temp[ni][nj]
						surround.append((ni, nj))
				cheat /= hcount
				for x, y in surround:
					temp[x][y] = cheat

		for i in range(ROWS):
			temp[i][0] = 22 + 273.15

		for i in range(ROWS):
			for j in range(COLS):
				if temp[i][j] > maxtemp[i][j]:
					failure = True
					break
			if failure:
				break
		if failure:
			print(temp)
			break
	if failure:
		failures += 1

print(failures)