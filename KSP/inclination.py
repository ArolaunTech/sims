import numpy as np
import matplotlib.pyplot as plt

changes = np.linspace(0, 90, num=1000)
I = changes * np.pi/180

x = 1 - np.cos(I)
frac = np.where(
	x < 2/9, 
	np.sqrt(2*x), 
	np.where(
		x < 1/2,
		4*np.sqrt(np.sqrt(2*x)-x)-2,
		2*(np.sqrt(2)-1)
	)
)

plt.ylim(0,1)
plt.xlim(0,90)
plt.gca().set_aspect(180/np.pi)
plt.title("Minimum Δv needed to change inclination from circular orbit")
plt.xlabel("Inclination change (degrees)")
plt.ylabel("Δv required (fraction of orbital velocity in circular orbit)")
plt.gca().fill_between(
	[0, np.arccos(7/9)*180/np.pi],
	0,
	1,
	color=(1,0.3,0,0.2)
)
plt.gca().fill_between(
	[np.arccos(7/9)*180/np.pi, 60],
	0,
	1,
	color=(0.15,0.6,0.3,0.2)
)
plt.gca().fill_between(
	[60, 90],
	0,
	1,
	color=(0.15,0.3,0.8,0.2)
)
plt.vlines([np.arccos(7/9)*180/np.pi, 60], 0, 1, colors="black", linestyles="dashed")
plt.text(
	np.arccos(7/9)*90/np.pi,
	0.7,
	"0° < θ < 38.94°\nDo inclination change\nin circular orbit",
	color="black",
	ha="center",
	fontsize="medium"
)
plt.text(
	np.arccos(7/9)*90/np.pi+30,
	0.15,
	"38.94° < θ < 60°\nRaise apoapsis\nand do inclination\nchange there",
	color="black",
	ha="center",
	fontsize="medium"
)
plt.text(
	75,
	0.2,
	"60° < θ < 90°\nBi-parabolic transfer",
	color="black",
	ha="center",
	fontsize="medium"
)
plt.plot(changes, frac, color="black")
plt.show()