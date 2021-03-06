
'''Script containing parameters and functions for Eady model scripts, following Vallis 2006.'''

from __future__ import division

import numpy as np
import cmath

################################# PHYSICAL PARAMETERS #############################################################

#Oceanic parameters

f0 = -1e-4 #(s^-1) f-plane approximation (southern hemisphere)
N = 2e-3 #(s^-1) buoyancy frequency (for SO, assuming uniform stratification)
H = 1e3 #(m) height of upper boundary
U_0 = 0.3 #(ms^(-1)) mean flow zonal velocity magnitude at upper boundary, faster to better represent eddy velocities in ACC
Lambda = U_0/H #(s^-1) uniform velocity shear in vertical direction
L = np.abs((N*H)/f0) #(m) typical length scale given by deformation radius
T = L/U_0 #(s) Eady timescale (about a day)
k = 1./L  #(m^-1) zonal wavenumber
l = (np.pi)/(2*L) #(m^-1) meridional wavenumber, defined for zonal channel
mu = L*np.sqrt(k**2 + l**2) #dimensionless parameter governing vertical wave struture

#Phase speed (ms^(-1)) (taking positive root for growth), Vallis 6.86
c = U_0/2. + (U_0/mu)*cmath.sqrt((mu/2. - (1./np.tanh(mu/2.)))*(mu/2. - np.tanh(mu/2.)))

sigma = k*c.imag #(s^-1) growth rate

################################# FUNCTIONS #############################################################

#Express all functions in terms of nondimensional variables x, y, z, t
#Dimensional variables are x^* = Lx, y^* = Ly, z^* = Hz and t^* = Tt
#Velocity perturbations are nondimensionalised by factors of (T/L), (or (T/H) for vertical).


#Real part of the vertical structure of the streamfunction perturbation, Vallis 6.93
def phi_r(z):
    return np.cosh(mu*z) - ((U_0*c.real)/(mu*(np.abs(c)**2)))*np.sinh(mu*z)

#Imaginary part of the vertical structure of the streamfunction perturbation, Vallis 6.93
def phi_i(z):
    return ((U_0*c.imag)/(mu*(np.abs(c)**2)))*np.sinh(mu*z)

#Derivative of phi_r(z) with respect to z
def dphi_rdz(z):
    return (mu/H)*np.sinh(mu*z) - ((U_0*c.real)/(H*(np.abs(c)**2)))*np.cosh(mu*z)

#Derivative of phi_i(z) with respect to z
def dphi_idz(z):
    return ((U_0*c.imag)/(H*(np.abs(c)**2)))*np.cosh(mu*z)

#Streamfunction perturbation, Vallis 6.77
def psiprime(x,y,z,t):
    return (phi_r(z)*np.cos(k*(L*x - c.real*T*t)) - phi_i(z)*np.sin(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t)

#Nondimensional zonal flow of basic state
def umeanflow(z):
    return (T/L)*Lambda*H*z


#Velocity perturbation amplitudes currently scaled by 100
#Note that as vprime is used in the wprime expression, vprime should be decreased accordingly so that
#the amplitude of wprime is scaled correctly

#Define nondimensional velocity perturbations (see notebook for calculations)
def uprime(x,y,z,t):
    return 100*(T/L)*l*(phi_r(z)*np.cos(k*(L*x - c.real*T*t)) - phi_i(z)*np.sin(k*(L*x - c.real*T*t)))*np.sin(l*L*y)*np.exp(sigma*T*t)

def vprime(x,y,z,t):
    return -100*(T/L)*k*(phi_r(z)*np.sin(k*(L*x - c.real*T*t)) + phi_i(z)*np.cos(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t)

def wprime(x,y,z,t):
    return -100*(T/H)*(f0/(N**2))*(((k*(c.real - Lambda*H*z)*dphi_rdz(z) - sigma*dphi_idz(z))*np.sin(k*(L*x - c.real*T*t)) + (k*(c.real - Lambda*H*z)*dphi_idz(z) + sigma*dphi_rdz(z))*np.cos(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t) - Lambda*(L/T)*0.01*vprime(x,y,z,t))
