'''Script containing parameters and functions for Eady model scripts, following Vallis 2006.'''

from __future__ import division

import numpy as np
import cmath

################################# PHYSICAL PARAMETERS #############################################################

#Oceanic parameters
f0 = -1e-4 #(s^-1) f-plane approximation (southern hemisphere)
N = 2e-3 #(s^-1) (King et al. 2012, JGR) buoyancy frequency (for SO, assuming uniform stratification)
#H = 1e3 #(m) height of upper boundary, meant to represent typical vertical extent of an oceanic eddy
H = 3e3
U_0 = 0.3 #(ms^(-1)) REFERENCE? mean flow zonal velocity magnitude at upper boundary
#U_0 = 0.1 #(ms^(-1)) considering difference in mean zonal flow from MITgcm simulations in first 1km depth
Lambda = U_0/H #(s^-1) uniform velocity shear in vertical direction
L = np.abs((N*H)/f0) #(m) typical length scale given by deformation radius
T = L/U_0 #(s) Eady timescale (about 2 days)
#k = 1./L  #(m^-1) zonal wavenumber
#l = (np.pi)/(2*L) #(m^-1) meridional wavenumber, defined for zonal channel
l = (np.pi)/(8*L)

alpha = 2e-4 #(K^-1) REFERENCE? thermal expansion coefficient
g = 9.8 #(ms^(-2)) gravitational acceleration
rho_0 = 1000 #(kg/m^3) REFERENCE? density of seawater
c_p = 4000 #(J/kg/K) REFERENCE? specific heat capacity of seawater

#Information to plot background theta distribution
theta0 = 278. #Reference temperature value, about 5 degrees celsius
dthetadz = (N**2)/(g*alpha) #Gill 6.17.24, buoyancy frequency definition
dthetady = (-f0*Lambda)/(g*alpha) #Gill 7.7.10, thermal wind

#Choose k based on chosen value of mu
def k_function(muvar):
    return np.sqrt(((muvar**2)/(L*L)) - l*l)

#Phase speed (ms^(-1)) (taking positive root for growth), Vallis 6.86
def c_function(muvar):
    return U_0/2. + (U_0/muvar)*cmath.sqrt((muvar/2. - (1./np.tanh(muvar/2.)))*(muvar/2. - np.tanh(muvar/2.)))

 #growth rate (s^-1)
def sigma_function(muvar):
    return k_function(muvar)*c_function(muvar).imag

################################# FUNCTIONS #############################################################

#Real part of the vertical structure of the streamfunction perturbation, Vallis 6.93
def phi_r_function(z,muvar):
    return np.cosh(muvar*z) - ((U_0*c_function(muvar).real)/(muvar*(np.abs(c_function(muvar))**2)))*np.sinh(muvar*z)

#Imaginary part of the vertical structure of the streamfunction perturbation, Vallis 6.93
def phi_i_function(z,muvar):
    return ((U_0*c_function(muvar).imag)/(muvar*(np.abs(c_function(muvar))**2)))*np.sinh(muvar*z)

#Derivative of phi_r(z) with respect to z
def dphi_rdz_function(z, muvar):
    return (muvar/H)*np.sinh(muvar*z) - ((U_0*c_function(muvar).real)/(H*(np.abs(c_function(muvar))**2)))*np.cosh(muvar*z)

#Derivative of phi_i(z) with respect to z
def dphi_idz_function(z, muvar):
    return ((U_0*c_function(muvar).imag)/(H*(np.abs(c_function(muvar))**2)))*np.cosh(muvar*z)

#Derivative of streamfunction with respect to z
def dpsi_dz_function(x,y,z,t,muvar):
    return (dphi_rdz_function(z,muvar)*np.cos(k_function(muvar)*(L*x - c_function(muvar).real*T*t)) - dphi_idz_function(z,muvar)*np.sin(k_function(muvar)*(L*x - c_function(muvar).real*T*t)))*np.cos(l*L*y)*np.exp(sigma_function(muvar)*T*t)


def vprime_function(x,y,z,t,muvar):
    return -k_function(muvar)*(phi_r_function(z,muvar)*np.sin(k_function(muvar)*(L*x - c_function(muvar).real*T*t)) + phi_i_function(z,muvar)*np.cos(k_function(muvar)*(L*x - c_function(muvar).real*T*t)))*np.cos(l*L*y)*np.exp(sigma_function(muvar)*T*t)

def wprime_function(x,y,z,t,muvar):
    return -(f0/(N**2))*(((k_function(muvar)*(c_function(muvar).real - Lambda*H*z)*dphi_rdz_function(z,muvar) - sigma_function(muvar)*dphi_idz_function(z,muvar))*np.sin(k_function(muvar)*(L*x - c_function(muvar).real*T*t)) + (k_function(muvar)*(c_function(muvar).real - Lambda*H*z)*dphi_idz_function(z,muvar) + sigma_function(muvar)*dphi_rdz_function(z,muvar))*np.cos(k_function(muvar)*(L*x - c_function(muvar).real*T*t)))*np.cos(l*L*y)*np.exp(sigma_function(muvar)*T*t) - Lambda*vprime_function(x,y,z,t,muvar))


#Conversion of APE density to KE density, Vallis 5.6.1
def APEconversion(x,y,z,t,muvar):
    return f0*wprime_function(x,y,z,t,muvar)*dpsi_dz_function(x,y,z,t,muvar)
