import numpy as np
from numpy.lib.scimath import sqrt
import mpmath
import plasmapy.mathematics as plmath
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
import sys
import fnmatch
import os

#read phi2 data from gkyl run
import postgkyl as pg

k=1.0
kperpsC = np.logspace(-2,0)

for beta in [0.1, 1.0, 10.0]:
    freqsC=[]

    for kperp in kperpsC:
        root = mpmath.findroot(lambda w: 
            beta*complex(w)**2/k**2*(1+w/np.sqrt(2)/k*plmath.plasma_dispersion_func(complex(w)/np.sqrt(2)/k))
            -(1+kperp**2+complex(w)/np.sqrt(2)/k*plmath.plasma_dispersion_func(complex(w)/np.sqrt(2)/k)),1-1j)
        freqsC.append(np.float64(root.real))

    plt.plot(kperpsC,freqsC,label=r"$\hat{\beta}=$"+str(beta))

    kperps = []
    freqs = []
    dirname = "beta-" + str(beta)
    for file in os.listdir(dirname):
        if fnmatch.fnmatch(file, "*_phi2_1.bp"):
            filename = file[:-4]
    
            data = pg.data.GData(dirname + "/" + filename)
            kperps.append(np.float(filename.split("_")[0].split("-")[1]))
            
            #put timesteps and phi2 values into lists
            times = data.peakGrid()[0]
            phi2s = data.peakValues()
            
            peaktimes = times[argrelextrema(phi2s, np.greater)[0]]
            diff = 0
            count = 0
            for i in np.arange(peaktimes.size-3):
                    diff = diff + peaktimes[i+1]-peaktimes[i]
                    count = count + 1
            period = 2*diff/count
            freqs.append(2*np.pi/period)
    
    kperps = np.array(kperps)
    freqs = np.array(freqs)
    # sort
    inds = kperps.argsort()
    kperps = kperps[inds]
    freqs = freqs[inds]
    
    if beta==10.0:
        plt.plot(kperps, freqs, "ok", label="Gkeyll")
    else:
        plt.plot(kperps, freqs, "ok")

plt.xscale("log")
plt.ylim(0,4)
plt.legend(loc="upper right")
plt.xlabel(r"$k_\perp \rho_s$")
plt.ylabel(r"$\omega/k_\parallel v_{te}$")

plt.savefig("alfven-kperp-scans.png")
plt.show()
