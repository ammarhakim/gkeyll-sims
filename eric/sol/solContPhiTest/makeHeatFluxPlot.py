from pylab import *
import numpy
import tables

from matplotlib import rcParams
import matplotlib.pyplot as plt

filenumEnd = 5

pathBase = '/Users/eshi/Research/gkeyllall/gkeyll-sims/eric/sol/solContPhiTest/'
fileBase = 'nocont_'
dataname = 'heatFluxAtEdge'

# Make a plot
fig = plt.figure(1, figsize=(9,4.75))

for idx in range(filenumEnd+1):
  filename = pathBase + fileBase + dataname + '_' + str(idx) + '.h5'
  fh = tables.openFile(filename)

  timeVals = fh.root.DataStruct.timeMesh.read()
  dataVals = fh.root.DataStruct.data.read()

  print dataVals.shape
  #plt.semilogx(timeVals, dataVals)

#plt.show()
