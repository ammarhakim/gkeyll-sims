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

kPerpSimList = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
tSimList = [2.980e-7, 1.327e-6, 2.061e-6, 2.549e-6, 2.652e-6, 2.607e-6, 2.544e-6, 2.479e-6, 2.376e-6, 2.280e-6,
    2.188e-6, 2.098e-6]
omegaSimList = zeros(len(kPerpSimList))

for index, kPerpRho in enumerate(kPerpSimList):
  omegaSimList[index] = 0.5*(2*pi/tSimList[index])/(kPar*vA)

plt.plot(kPerpRhoList, exactFreqList,'r-',label='Exact')
plt.plot(kPerpSimList, omegaSimList,'b-o',label='Sim')

#plt.xlim(kPerpRhoList[0],kPerpRhoList[-1])
plt.ylim(0,8)
plt.xlabel(r'$k_\perp \rho_s$')
plt.ylabel(r'$\omega/(k_\parallel v_A)$')
plt.legend(loc='upper right')
#plt.autoscale(enable=True,axis='y',tight=True)
plt.title(r'Scan at 16X32V, $\beta_e$ = 0.01')
plt.savefig('kaw_kPerpScan.pdf')
plt.show()
#plt.close()
