"""#===== Lattice Boltzmann Method =====

Uses the BGK approximation in order to collide particles.

All units are SI in order to avoid unit confusion.

"""#=====            END           =====

#Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.colors as colors
from timemodule import timer

a = timer()
a.tagtime()

#Simulation constants
ids = np.arange(9)
weights = np.array([1/36, 1/9, 1/36, 1/9, 4/9, 1/9, 1/36, 1/9, 1/36])
velr = np.array([-1, -1, -1, 0, 0, 0, 1, 1, 1])
velc = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1])

#Simulation settings
s = 0.01 #Side length of square

rows = 100
cols = 100

R, C = np.meshgrid(range(rows), range(cols))

dt = 0.1 #Timestep
tEnd = 100

tau = 0.4 #Collision timescale

N = int(tEnd / dt) #int() is used here instead of the // floor division operator.

#Initial conditions
F = np.full((cols, rows, 9), [0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)

rho = np.sum(F,2)
vr = np.where(rho==0, 0, np.sum(F*velr, 2) / rho)
vc = np.where(rho==0, 0, np.sum(F*velc, 2) / rho)

block = R > 100

fanr = np.logical_and(np.logical_and(R > 30, R < 70), C == 0)

history = np.array([np.where(block, np.nan, rho).T], dtype=np.float64)
histvelr = np.array([np.where(block, np.nan, vr).T], dtype=np.float64)
histvelc = np.array([np.where(block, np.nan, vc).T], dtype=np.float64)

print(a.elapsed())
a.tagtime()

#Simulation main loop
for it in range(N):
	#Stream
	for i in range(9):
		F[:,:,i] = np.roll(F[:,:,i], (velr[i], velc[i]), axis=(1,0))

	#Set boundaries
	bdry = F[block, :][:, [8,7,6,5,4,3,2,1,0]]

	#Calculate variables
	rho = np.sum(F, 2)
	vr = np.where(rho==0, 0, np.sum(F*velr, 2) / rho)
	vc = np.where(rho==0, 0, np.sum(F*velc, 2) / rho)

	#vr = 0.005*(50 - R)
	#vc = 0.005*(50 - C)

	np.clip(vr, -1, 1, out=vr)
	np.clip(vc, -1, 1, out=vc)
	"""
	dist = vr*vr + vc*vc
	np.sqrt(dist, out=dist)
	np.clip(dist, 1/np.sqrt(3), None, dist)
	dist *= np.sqrt(3)

	vr /= dist
	vc /= dist
	"""

	#Collide
	eq = np.zeros(F.shape)
	for i, cr, cc, w in zip(ids, velr, velc, weights):
		eq[:,:,i] = rho * w * (1 + 3 * (cr*vr+cc*vc) + 4.5 * (cr*vr+cc*vc)**2 - 1.5 * (vr*vr+vc*vc)**2)

	s = np.sum(eq, 2)

	for i in range(9):
		eq[:,:,i] *= np.where(s==0, 0, rho / s)

	F += tau * (eq - F)

	F[block, :] = bdry

	#Manipulate
	bdry = F[fanr, :]
	F[fanr, :] = 0.0
	F[fanr, 5] = np.sum(bdry, 1)

	history = np.append(history, np.array([np.where(block, np.nan, rho).T]), axis=0)
	histvelr = np.append(histvelr, np.array([np.where(block, np.nan, vr).T]), axis=0)
	histvelc = np.append(histvelc, np.array([np.where(block, np.nan, vc).T]), axis=0)

print(a.elapsed())
a.tagtime()

RV = np.arange(0, rows, 10)
CV = np.arange(0, cols, 10)
RV2, CV2 = np.meshgrid(CV, RV)

history = history.tolist()

curr = cm.get_cmap("inferno").copy()
curr = colors.ListedColormap(curr(np.linspace(0.1,1,230)))
curr.set_bad(color="black")

div = cm.get_cmap("bwr").copy()
div.set_bad(color="black")

fig, ax = plt.subplots()
im = plt.imshow(history[0], cmap=curr, vmin=0, vmax=1, animated=True)
plt.title("Fluid density")
plt.colorbar()
qax = plt.quiver(RV2, CV2, histvelc[0][RV][:,CV].tolist(), histvelr[0][RV][:,CV].tolist(), units="xy", angles="xy", scale=0.01, scale_units="xy")

def animate(i):
	im.set_array(history[i])
	qax.set_UVC(histvelc[i][RV][:,CV].tolist(), histvelr[i][RV][:,CV].tolist())
	return im,

ani = animation.FuncAnimation(fig, animate, interval=50, blit=True, frames=N)
ani.save('simulation07112022.gif', writer='imagemagick', fps=20)
print(a.elapsed())
plt.show()