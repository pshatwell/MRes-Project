import numpy as np
from sympy import *

################################# LIST OF PHYSICAL PARAMETERS #############################################################

#Typical parameters for ocean, for most unstable (fastest growing) Eady wave

L = 1e4 #(m) typical length scale
f0 = -1e-4 #(s^-1) f-plane approx. for Coriolis parameter (southern hemisphere)
T = 10*np.abs(1./f0) #(s) typical time scale
N = 1e-3 #(s^-1) buoyancy frequency
k = 1e-4 #(m^-1) zonal wavenumber of growing wave
H_R = np.abs(f0/(N*k)) #(m) Rossby height
H = 0.8031*H_R #(m) positions of boundaries at H and -H for most unstable wave
Hratio = np.abs(H/H_R)
U = 0.1 #(ms^(-1)) mean flow zonal velocity magnitude at boundaries H and -H
shear = U/H #(s^-1) velocity shear
sigma_max = np.abs(0.3098*(f0/N)*shear) #maximum growth rate
alpha = 2e-4 #(K^-1) thermal expansion coefficient
g = 9.8 #(ms^(-2)) gravitational acceleration

########################## DEFINE STREAMFUNCTION AND DERIVATIVES ####################################################################

#def phiprime(x,z,t): #Define function for streamfunction perturbation
#    return (np.cos(k*x)*np.divide(np.sinh(z/H_R), np.sinh(Hratio)) + np.sin(k*x)*np.divide(np.cosh(z/H_R), np.cosh(Hratio)))*np.exp(sigma_max*t)

x,z,t = symbols('x z t')

#Define streamfunction perturbation for most unstable solution, Gill 13.3.15
phiprime = (cos(k*L*x)*(sinh(z*Hratio)/sinh(Hratio)) + sin(k*L*x)*(cosh(z*Hratio)/cosh(Hratio)))*exp(sigma_max*t)

#use Sympy to find expressions for the partial derivatives
phiprime_x = (1./L)*diff(phiprime, x)
phiprime_z = (1./H)*diff(phiprime, z)
phiprime_t = np.abs(f0)*diff(phiprime, t)

phiprime_zx = diff(phiprime_x, z)
phiprime_zt = diff(phiprime_t, z)
phiprime_xx = diff(phiprime_x, x)
