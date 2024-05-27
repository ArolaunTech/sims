import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from noise import snoise2
from scipy.signal import convolve2d

nf = np.vectorize(snoise2)

path = 'camera.png'
data = np.float64(np.array(Image.open(path)))/255
if len(data.shape) == 3:
	data = data[:,:,0]

x = np.arange(0, 1001)
X, Y = np.meshgrid(x,x)
X -= 500
Y -= 500
data = np.zeros(X.shape)
#data = r
data = 200 - np.sqrt(X**2 + Y**2)
data += 10 * np.random.random(X.shape)
data += 5 * nf(X/100, Y/100)
data += 2 * nf(X/50, Y/50)
data += 1 * nf(X/20, Y/20)
#data = np.where(X**2 + Y**2 <= 10000* r**2, 1, 0)

data *= 0.8+0.4*nf(X/100,Y/100+5)

data = np.maximum(np.zeros(data.shape), data)
data = np.minimum(np.ones(data.shape), data)

xcircle = np.arange(-5, 6)
Xcircle, Ycircle = np.meshgrid(xcircle, xcircle)
circle = np.where(Xcircle**2 + Ycircle**2 <= 20, 0.05, 0)

ft = np.fft.ifftshift(data)
ft = np.fft.fft2(ft)
ft = np.fft.fftshift(ft)
ft = abs(ft)

ft = np.maximum(np.zeros(data.shape), ft)
ft /= np.max(ft)

ft = convolve2d(ft, circle, mode='same')
ft *= 10

fig, ax = plt.subplots(1,2)
ax[0].imshow(data, cmap='inferno')
ax[1].imshow(np.dstack((ft, ft, 0.7*ft)))
plt.show()

imgout = Image.fromarray(np.uint8(255*np.minimum(1,np.dstack((ft, ft, 0.7*ft)))))
imgout.save('kerbol.png')
