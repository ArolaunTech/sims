import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import matplotlib.cm as cm

def xval(time, fxc1, fxc2, fxc3):
	return fxc2*np.cos(time)-fxc1+fxc3*np.sin(time)

delta = 1000
angles = 2*np.pi*np.arange(0,1,1/delta)

yeartimes = np.arange(0,2*np.pi,2*np.pi/365)

lats = [-70,-66.5,-65,-60,-38,0,38,60,65,66.5,70]

def calc_line(latitude):
	fxc1 = np.sin(latitude)*np.sin(np.pi/180 * 23.5)
	fxc2 = np.cos(latitude)*np.cos(np.pi/180 * 23.5)
	fxc3 = np.cos(latitude)
	out = []
	for yeartime in yeartimes:
		out.append(np.sum(xval(angles, fxc1*np.cos(yeartime), fxc2*np.cos(yeartime), fxc3*np.sin(yeartime))>0)*0.024)
	return out

solstice = dt.datetime.strptime('12/21','%m/%d').date()
dates = [solstice+dt.timedelta(days=i) for i in range(365)]

fig, ax = plt.subplots()

plt.title("Day lengths for different latitudes and times")
plt.xlabel('Month')
plt.ylabel('Day Length (Hours)')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

ax.set_yticks(np.linspace(0,25,26), minor=True)
ax.set_yticks(np.linspace(0,25,6))

ax.grid(which='both')

# Or if you want different settings for the grids:
ax.grid(which='minor', alpha=0.5)
ax.grid(which='major', alpha=1)

np.random.seed(431085341)

def sigmoid(x):
	return 1/(1+np.exp(-x))

cmap = cm.get_cmap('turbo')

for i, lat in enumerate(lats):
	lab = str(abs(lat)) + '°'
	if lat >= 0:
		lab += 'N'
	else:
		lab += 'S'
	if lat == 66.5:
		lab = 'Arctic Circle'
	if lat == -66.5:
		lab = 'Antarctic Circle'
	col = cmap(i/10)
	plt.plot(dates, calc_line(np.pi/180*lat), color=col, label=lab)

plt.legend(loc='upper right')

plt.show()