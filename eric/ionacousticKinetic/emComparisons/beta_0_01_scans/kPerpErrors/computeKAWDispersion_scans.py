from numpy import *
from scipy import special
from scipy import optimize
import matplotlib.pyplot as plt
import math
from matplotlib import rcParams
from matplotlib import rc

rcParams['legend.fontsize'] = 22
font = {'family' : 'sans-serif',
        'sans-serif' : ['Arial'],
        'size'   : 22}
rc('font', **font)

def plasmaDisp(z):
    return 1j*sqrt(pi)*exp(-z**2)*(1+special.erf(1j*z))

def derivPlasmaDisp(z):
    return -2*(1+z*plasmaDisp(z))

def eps(z,beta_e_val):
    return (m_i/m_e*beta_e_val*z**2-1)*(1+z*plasmaDisp(z)) - kPerpRho**2

def derivEps(z,beta_e_val):
    return (2*m_i/m_e*beta_e_val*z)*(1+z*plasmaDisp(z)) + (m_i/m_e*beta_e_val*z**2-1)*(plasmaDisp(z)+z*derivPlasmaDisp(z))

# Number of points to calculate damping rate at
nPoints = 100;
kPerpRhoList = linspace(0.01, 1.0, nPoints)
exactFreqList = zeros(nPoints)

m_i = 1.672621777e-27
m_e = 9.10938188e-31
mu0 = 4e-7*pi
kPar = 0.5
B = 1
Te0 = 250
eV = 1.602176565e-19
nSim = 9.947e19
beta_e = 2*mu0*Te0*eV*nSim/B**2
vA = B/sqrt(mu0*m_i*nSim)

tol = 1e-8
for index, kPerpRho in enumerate(kPerpRhoList):
  # Initial guess for z0 = omega/(k*sqrt(2)*vTe) using approximate expression for wave frequency
  z0 = vA/(sqrt(2*Te0*eV/m_e)*sqrt(1+2/beta_e*m_e/m_i*kPerpRho**2))

  z0 = optimize.newton(eps,z0,derivEps,(beta_e,),tol,10000)
  exactFreqList[index] = fabs(z0.real*sqrt(2*Te0*eV/m_e)/vA);

kPerpSimList = [0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
tSimList = [2.802e-6, 2.861e-6, 2.867e-6, 2.861e-6, 2.829e-6, 2.759e-6, 2.681e-6, 2.589e-6, 2.492e-6, 2.392e-6, 2.294e-6,
    2.199e-6, 2.106e-6]
tSimListErrors = [2.980e-7, 8.576e-7, 1.327e-6, 2.061e-6, 2.549e-6, 2.652e-6, 2.607e-6, 2.544e-6, 2.479e-6, 2.376e-6, 2.280e-6,
    2.188e-6, 2.098e-6]
omegaSimList = zeros(len(kPerpSimList))
omegaSimListErrors = zeros(len(kPerpSimList))

for index, kPerpRho in enumerate(kPerpSimList):
  omegaSimList[index] = 0.5*(2*pi/tSimList[index])/(kPar*vA)
  omegaSimListErrors[index] = 0.5*(2*pi/tSimListErrors[index])/(kPar*vA)

# Make a plot
fig = plt.figure(1, figsize=(10,8))

plt.semilogx(kPerpRhoList, exactFreqList,'r-',label='Exact')
plt.semilogx(kPerpSimList, omegaSimList,'bo',label='Sim')
plt.semilogx(kPerpSimList, omegaSimListErrors,'g-o',label='0.1% Error in Density')

#plt.xlim(kPerpRhoList[0],kPerpRhoList[-1])
plt.ylim(0,10)
plt.xlabel(r'$k_\perp \rho_s$')
plt.ylabel(r'$\omega/(k_\parallel v_A)$')
plt.legend(loc='upper right')
ax = fig.gca()
ax.set_xscale('log', basex=10, subsx=arange(2,9,1))
ax.xaxis.grid(True, which='both')
#plt.autoscale(enable=True,axis='y',tight=True)
#plt.title(r'Scan at 16X32V, $\beta_e$ = 0.01')
plt.savefig('kaw_kPerpScan.pdf',bbox_inches='tight', \
                        pad_inches=0.1)
plt.show()
#plt.close()
