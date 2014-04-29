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

pathBase = '/Users/eshi/Research/gkeyllall/gkeyll-sims/eric/sol/solElectrostaticPaper/'
fileBase = 'es1_'
dataname = 'heatFluxAtEdge'

# Make a plot
fig = plt.figure(1, figsize=(10,8))

for idx in range(1,filenumEnd+1):
  filename = pathBase + fileBase + dataname + '_' + str(idx) + '.h5'
  print filename
  fh = tables.openFile(filename)

  timeVals = fh.root.DataStruct.timeMesh.read()*1e6
  dataVals = fh.root.DataStruct.data.read()

  plt.semilogx(timeVals, dataVals[:,0],'b')
  plt.semilogx(timeVals, dataVals[:,1],'r')
  plt.semilogx(timeVals, dataVals[:,2],'k')

plt.ylim(0, 6e9)
plt.xlim(0.5, 400)
plt.legend(['Total','Ions','Electrons'],loc='upper left')
plt.xlabel(r'time ($\mu$s)')
plt.ylabel(r'Q (W/m$^2$)')
plt.grid()
ax = fig.gca()
ax.set_xscale('log', basex=10, subsx=np.arange(2,9,1))
ax.xaxis.grid(True, which='both')

plt.savefig('kinetic-elc-heat-flux1.pdf', bbox_inches='tight')
plt.show()
