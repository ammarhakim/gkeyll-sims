import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import sys

#read phi2 data from gkyl run
import postgkyl as pg
filename = sys.argv[1]
data = pg.data.GData(filename)
beta = np.float(filename.split("_")[0].split("-")[1])

expected = 1.0/np.sqrt(beta+.01)

#put timesteps and phi2 values into lists
times = data.peakGrid()[0]
phi2s = data.peakValues()

peaktimes = times[argrelextrema(phi2s, np.greater)[0]]
peaktimes = np.insert(peaktimes, 0, 0)
diff = 0
count = 0
for i in np.arange(peaktimes.size-1):
        diff = diff + peaktimes[i+1]-peaktimes[i]
        count = count + 1
period = 2*diff/count
freq = 2*np.pi/period
print("%f\t%f\t%f" %( beta, freq, expected))
