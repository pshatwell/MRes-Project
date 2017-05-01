

'''Script containing parameters and functions for Eady model scripts, following Vallis 2006.'''


import numpy as np
import cmath


################################# PHYSICAL PARAMETERS #############################################################

#Oceanic parameters, Vallis p. 269

f0 = -1e-4 #(s^-1) f-plane approximation (southern hemisphere)
N = 1e-2 #(s^-1) buoyancy frequency (assuming uniform stratification)
H = 1e3 #(m) height of upper boundary
U_0 = 0.1 #(ms^(-1)) mean flow zonal velocity magnitude at upper boundary
Lambda = U_0/H #(s^-1) uniform velocity shear in vertical direction
L = np.abs((N*H)/f0) #(m) typical length scale given by deformation radius
T = L/U_0 #(s) Eady timescale
k = 1e-5 #(m^-1) zonal wavenumber
#l = (np.pi)/L #(m^-1) meridional wavenumber, defined for zonal channel
l = 1e-10
mu = L*np.sqrt(k**2 + l**2) #dimensionless parameter governing vertical wave struture

alpha = 2e-4 #(K^-1) thermal expansion coefficient
g = 9.8 #(ms^(-2)) gravitational acceleration

#Phase speed (taking positive root - IS THIS JUSTIFIED??), Vallis 6.86
c = U_0/2. + (U_0/mu)*cmath.sqrt((mu/2. - (1./np.tanh(mu/2.)))*(mu/2. - np.tanh(mu/2.)))

sigma = k*c.imag #(s^-1) growth rate


################################# FUNCTIONS #############################################################

#Express all functions in terms of nondimensional variables x, y, z, t
#Dimensional variables are x^* = Lx, y^* = Ly, z^* = Hz and t^* = Tt


#Real part of the vertical structure of the streamfunction perturbation, Vallis 6.93
def phi_r(z):
    return np.cosh(mu*z) - ((U_0*c.real)/(mu*(np.abs(c)**2)))*np.sinh(mu*z)

#Derivative of phi_r(z) with respect to z
def dphi_rdz(z):
    return (mu/H)*np.sinh(mu*z) - ((U_0*c.real)/(H*(np.abs(c)**2)))*np.cosh(mu*z)

#Streamfunction perturbation, Vallis 6.77
def psiprime(x,y,z,t):
    return phi_r(z)*np.sin(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(sigma*T*t)

#Nondimensional zonal flow of basic state
def umeanflow(z):
    return (T/L)*Lambda*H*z


#Define nondimensional velocity perturbations (see notebook for calculations)
def uprime(x,y,z,t):
    return -(T/L)*l*phi_r(z)*np.cos(l*L*y)*np.cos(k*(L*x - c.real*T*t))*np.exp(sigma*T*t)

def vprime(x,y,z,t):
    return -(T/L)*k*phi_r(z)*np.sin(l*L*y)*np.sin(k*(L*x - c.real*T*t))*np.exp(sigma*T*t)

#SIGN OF f0 HERE??
def wprime(x,y,z,t):
    #return (T/H)*(f0/(N**2))*dphi_rdz(z)*np.sin(l*L*y)*((Lambda*H*z - c.real)*k*np.sin(k*(L*x - c.real*T*t)) - sigma*np.cos(k*(L*x - c.real*T*t)))*np.exp(sigma*T*t)
    return (T/H)*(f0/(N**2))*dphi_rdz(z)*np.sin(l*L*y)*((Lambda*H*z - c.real)*k*np.sin(k*(L*x - c.real*T*t)) - sigma*np.cos(k*(L*x - c.real*T*t)) - Lambda*vprime(x,y,z,t))*np.exp(sigma*T*t)
