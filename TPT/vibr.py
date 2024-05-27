import matplotlib.pyplot as plt
import numpy as np

freq_array = [[0 for j in range(7)] for i in range(7)]
iters = 1e+6
for i in range(int(iters)):
	rng = np.random.randint(2**32) - 2**31
	rng = rng >> 4
	col = rng%7
	row = (rng>>3)%7
	freq_array[row][col] += 1/iters

plt.imshow(freq_array, vmin=0, cmap='inferno')
plt.colorbar()
plt.show()