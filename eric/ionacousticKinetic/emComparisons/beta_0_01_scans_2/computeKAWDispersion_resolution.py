# At one value of beta_e, plot resolution scan frequencies
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
m_i = 1.672621777e-27
m_e = 9.10938188e-31
mu0 = 4e-7*pi
kPar = 0.5
B = 1
Te0 = 250
eV = 1.602176565e-19

tol = 1e-8

kPerpRho = 0.01
nSim = 9.947e19

beta_e = 2*mu0*Te0*eV*nSim/B**2
vA = B/sqrt(mu0*m_i*nSim)
# Initial guess for z0 = omega/(k*sqrt(2)*vTe) using approximate expression for wave frequency
z0 = vA/(sqrt(2*Te0*eV/m_e)*sqrt(1+2/beta_e*m_e/m_i*kPerpRho**2))
# Root finder call
z0 = optimize.newton(eps,z0,derivEps,(beta_e,),tol,10000)
# Store exact result
exactFreq = fabs(z0.real*sqrt(2*Te0*eV/m_e)/vA);
# Lame way to create list of three elements of same value
exactFreqList = [exactFreq, exactFreq, exactFreq]

# Data measured from simulation
vResList = [16, 32, 64]
tSimList = [4.002e-7, 4.018e-7, 4.020e-7]
omegaSim = zeros(len(vResList))
# Compute normalized simulation frequencies
for index, tSim in enumerate(tSimList):
  # Really normalized to k*vA
  # Divide by 2 to compare with wave freq
  omegaSim[index] = 0.5*(2*pi/tSim)/(kPar*B/sqrt(mu0*nSim*m_i))

plt.plot(vResList, exactFreqList,'b-o',label='Exact')
plt.plot(vResList, omegaSim,'r-o',label='16X')

plt.xlabel('Total P-Space Cells')
plt.ylabel(r'$\omega/(k_\parallel v_A)$')
plt.legend(loc='right')
plt.title(r'Resolution Scan for $k_\perp \rho_s$=0.01, $\beta_e$=0.01')
plt.ylim([0,8])
#plt.autoscale(enable=True,axis='y',tight=True)
plt.savefig('kaw_0_01_scan.pdf')
plt.show()
#plt.close()
