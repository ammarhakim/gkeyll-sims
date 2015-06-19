from pylab import *
import numpy
import tables

from matplotlib import rcParams
from matplotlib import rc
import matplotlib.pyplot as plt

rcParams['lines.linewidth'] = 4
rcParams['xtick.major.pad'] ='8'
rcParams['ytick.major.pad'] ='8'
font = {'family' : 'sans-serif',
        'sans-serif' : ['Arial'],
        'size'   : 22}
rc('font', **font)

filenumEnd = 5

pathBase = '/Users/eshi/Research/gkeyllall/gkeyll-sims/eric/sol/solWithMu/pitts/'
fileBase = 'esWithCollisions_'
dataname = 'sheathCoefficients' # To get gamma, temperatures
dataname2 = 'cutoffV' # To get sheath potential
# Physical Parameters
elcMass = 9.10938215e-31
eV = 1.602176487e-19

# Make a plot
#fig = plt.figure(1, figsize=(10,8))
f, axarr = plt.subplots(3, sharex=True, figsize=(8,14))

for idx in range(1,filenumEnd+1):
  filename = pathBase + fileBase + dataname + '_' + str(idx) + '.h5'
  fh = tables.openFile(filename)

  timeVals = fh.root.DataStruct.timeMesh.read()*1e6
  dataVals = fh.root.DataStruct.data.read()

  axarr[0].semilogx(timeVals, dataVals[:,0],'b') # (Total)
  axarr[0].semilogx(timeVals, dataVals[:,1],'r') # (Ions)
  axarr[0].semilogx(timeVals, dataVals[:,2],'k') # (Electrons)
  axarr[1].semilogx(timeVals, dataVals[:,4],'r') # (Ion Temp)
  axarr[1].semilogx(timeVals, dataVals[:,3],'k') # (Electron Temp)
    
  filename = pathBase + fileBase + dataname2 + '_' + str(idx) + '.h5'
  fh = tables.openFile(filename)
  timeVals = fh.root.DataStruct.timeMesh.read()*1e6
  dataVals = fh.root.DataStruct.data.read()
  axarr[2].semilogx(timeVals, 0.5*elcMass/eV*dataVals[:,1]**2,'b') # (Sheath Potential)

axarr[0].legend(['Total','Ions','Electrons'],loc='upper left')
#ax1 = plt.subplot(211)
axarr[0].set_ylim(0,20)
axarr[1].set_ylim(0,1600)
axarr[1].set_xlim(0.5, 350)
axarr[0].set_title('Sheath Transmission Coefficient')
axarr[1].set_title('Temperature at Edge (eV)')
axarr[2].set_title('Sheath Potential (eV)')
#plt.xlabel(r'time ($\mu$s)')
#plt.ylabel(r'$\gamma$')
#plt.grid()
axarr[2].set_xscale('log', basex=10, subsx=np.arange(2,9,1))
axarr[0].xaxis.grid(True, which='both')
axarr[1].xaxis.grid(True, which='both')
axarr[2].xaxis.grid(True, which='both')

plt.savefig('sheath-coefficient.pdf', bbox_inches='tight')
plt.show()
