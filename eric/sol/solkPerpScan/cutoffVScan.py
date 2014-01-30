from pylab import *
import numpy
import tables

from matplotlib import rcParams
from matplotlib import rc
import matplotlib.pyplot as plt

rcParams['lines.linewidth']= 4
rcParams['xtick.major.pad']='8'
rcParams['ytick.major.pad']='8'
font = {'family' : 'sans-serif',
        'sans-serif' : ['Arial'],
        'size'   : 22}
rc('font', **font)

filenumEnd = 5

pathBase = '/Users/eshi/Research/gkeyllall/gkeyll-sims/eric/sol/solkPerpScan/'
fileList = ['0_02','0_1','0_2']
dataname = 'cutoffVelocities'

# Physical parameters
elcMass = 9.10938291e-31
eV = 1.60217657e-19

col = ['b', 'r', 'k', 'c', 'm']

# Make a plot
fig = plt.figure(1, figsize=(10,8))
plt.gca().set_color_cycle(['red','blue','green'])

for fileIndex in range(len(fileList)):
  fileBase = fileList[fileIndex]
  timeVals = []
  phiS = []
  
  for idx in range(1,filenumEnd+1):
    filename = pathBase + fileBase + '_' + dataname + '_' + str(idx) + '.h5'
    print filename
    fh = tables.openFile(filename)

    timeValsRead = fh.root.DataStruct.timeMesh.read()*1e6
    dataValsRead = fh.root.DataStruct.data.read()

    timeVals = append(timeVals, timeValsRead[:,0])
    phiS = append(phiS, 0.5*elcMass*dataValsRead[:,1]**2/eV)
  plt.semilogx(timeVals, phiS, col[fileIndex], label=fileBase.replace('_','.'))

#plt.ylim(0, 6e9)
plt.xlim(0.5, 400)
plt.legend(loc='lower center')
#plt.legend(['Total','Ions','Electrons'],loc='upper left')
plt.xlabel(r'time ($\mu$s)')
plt.ylabel(r'$\phi_s (eV)$')
plt.grid()
ax = fig.gca()
ax.set_xscale('log', basex=10, subsx=np.arange(2,9,1))
ax.xaxis.grid(True, which='both')
plt.savefig('sheathPotential.pdf', bbox_inches='tight')

plt.show()
