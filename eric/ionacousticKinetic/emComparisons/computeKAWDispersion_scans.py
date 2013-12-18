from numpy import *
from scipy import special
from scipy import optimize
import matplotlib.pyplot as plt
import math
from matplotlib import rc

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
beta_e_list = logspace(-4, -1, nPoints)
approxFreqList = zeros(nPoints);

m_i = 1.672621777e-27
m_e = 9.10938188e-31
mu0 = 4e-7*pi
kPar = 0.5
B = 1
Te0 = 250
eV = 1.602176565e-19

tol = 1e-4

kPerpRhoList = [0.1, 0.2, 0.6]
exactFreqList = []

for listIndex, kPerpRho in enumerate(kPerpRhoList):
  # Initial guess for z0 = omega/(k*sqrt(2)*vTe) using approximate expression for wave frequency
  n0 = B**2*beta_e_list[0]/(2*mu0*Te0*eV)
  vA = B/sqrt(mu0*m_i*n0)
  z0 = vA/(sqrt(2*Te0*eV/m_e)*sqrt(1+2/beta_e_list[0]*m_e/m_i*kPerpRho**2))
  
  freqList = zeros(nPoints)

  for index, beta_e in enumerate(beta_e_list):
      z0 = optimize.newton(eps,z0,derivEps,(beta_e,),tol,10000)

      nVal = B**2*beta_e/(2*mu0*Te0*eV)
      vA = B/sqrt(mu0*m_i*nVal)
      freqList[index] = fabs(z0.real*sqrt(2*Te0*eV/m_e)/vA);

      # Compute approximate solution using formula
      #approxFreqList[index] = 1/sqrt(1 + 2/beta_e_list[index]*m_e/m_i*kPerpRho**2)
  exactFreqList.append(freqList[:])

# Build list of simulation-derived parameters
nSim = [9.947e17, 9.947e18, 9.947e19, 9.947e20]
tSimList = []
# kPerpRho = 0.1
tSimList.append([2.962e-7, 8.851e-7, 2.336e-6, 3.683e-6])
# kPerpRho = 0.2
tSimList.append([3.206e-7, 8.970e-7, 2.657e-6, 5.972e-6])
# kPerpRho = 0.6
tSimList.append([4.382e-7, 8.697e-7, 2.495e-6, 7.403e-6])
omegaSimList = []
betaSimList = []

for listIndex, kPerpRho in enumerate(kPerpRhoList):
  omegaSim = zeros(len(nSim))
  betaSim = zeros(len(nSim))
  tSim = tSimList[listIndex]

  for index, n_val in enumerate(nSim):
    # Really normalized to k*vA
    # Divide by 2 to compare with wave freq
    omegaSim[index] = 0.5*(2*pi/tSim[index])/(kPar*B/sqrt(mu0*nSim[index]*m_i))
    betaSim[index] = nSim[index]*2*mu0*Te0*eV/(B**2)
  omegaSimList.append(omegaSim[:])
  betaSimList.append(betaSim[:])

plt.semilogx(beta_e_list, exactFreqList[0],'r-',label='Exact 0.1')
plt.semilogx(beta_e_list, exactFreqList[1],'b-',label='Exact 0.2')
plt.semilogx(beta_e_list, exactFreqList[2],'g-',label='Exact 0.6')
#plt.semilogx(beta_e_list, approxFreqList,'b-',label='Approx')
plt.semilogx(betaSimList[0], omegaSimList[0],'r-o',label='Sim 0.1')
plt.semilogx(betaSimList[1], omegaSimList[1],'b-o',label='Sim 0.2')
plt.semilogx(betaSimList[2], omegaSimList[2],'g-o',label='Sim 0.6')

plt.xlim(betaSimList[0][0],betaSimList[0][-1])
plt.xlabel(r'$\beta_e$')
plt.ylabel(r'$\omega/(k_\parallel v_A)$')
plt.legend(loc='upper left')
plt.autoscale(enable=True,axis='y',tight=True)
plt.savefig('kawKPerpScan.pdf')
plt.show()
#plt.close()
