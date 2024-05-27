import wave
import numpy as np
import matplotlib.pyplot as plt
import struct

#https://stackoverflow.com/questions/53460512/python-waveform-to-wav-file-converter
def signal_to_wav(signal, fname, Fs):
    """Convert a numpy array into a wav file.

     Args
     ----
     signal : 1-D numpy array
         An array containing the audio signal.
     fname : str
         Name of the audio file where the signal will be saved.
     Fs: int
        Sampling rate of the signal.

    """
    #print(*signal)
    data = struct.pack('<' + ('h'*len(signal)), *signal)
    wav_file = wave.open(fname, 'wb')
    wav_file.setnchannels(1)
    wav_file.setsampwidth(2)
    wav_file.setframerate(Fs)
    wav_file.writeframes(data)
    wav_file.close()

secondsLength = 2
sampleFreq = 40000
samples = secondsLength*sampleFreq
data = np.zeros(samples)
seconds = np.indices((samples,))[0]/sampleFreq

decay = 0.65
baseStart = 340
baseEnd = 320
baseMid = 335

c = (baseEnd - baseMid)/(baseMid - baseStart)
b = (baseMid - baseStart)/(c - 1)
a = baseStart - b

#base = a + b * c**(2*seconds/secondsLength)
freqSelect = 1-(4*seconds*seconds*np.exp(2-4*seconds))
soundDecay = np.exp(-seconds)
base = 262 + 10 * freqSelect
for i in range(1, 20):
    harmonic = np.sin(2*np.pi*seconds*base*i)
    harmonicMagnitude = np.power(decay, i-1)
    harmonicDecay = (soundDecay)**(2**(i-1))
    data += harmonicMagnitude*harmonic*harmonicDecay

data /= np.max(data)
data *= 25000
data = np.int64(data)

plt.plot(seconds, data)
plt.show()

signal_to_wav(data, 'hum.wav', sampleFreq)