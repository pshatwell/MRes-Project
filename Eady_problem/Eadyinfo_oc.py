
'''Script containing parameters and functions for Eady model scripts, following Vallis 2006.'''

from __future__ import division

import numpy as np
import cmath

################################# PHYSICAL PARAMETERS #############################################################

#Oceanic parameters
f0 = -1e-4 #(s^-1) f-plane approximation (southern hemisphere)
N = 2e-3 #(s^-1) (King et al. 2012, JGR) buoyancy frequency (for SO, assuming uniform stratification)
#H = 1e3 #(m) height of upper boundary, meant to represent typical vertical extent of an oceanic eddy
H = 2985.
U_0 = 0.3 #(ms^(-1)) REFERENCE? mean flow zonal velocity magnitude at upper boundary
#U_0 = 0.1 #(ms^(-1)) considering difference in mean zonal flow from MITgcm simulations in first 1km depth
Lambda = U_0/H #(s^-1) uniform velocity shear in vertical direction
L = np.abs((N*H)/f0) #(m) typical length scale given by deformation radius
T = L/U_0 #(s) Eady timescale (about 2 days)
#k = 1./L  #(m^-1) zonal wavenumber
#l = (np.pi)/(2*L) #(m^-1) meridional wavenumber, defined for zonal channel
l = (np.pi)/(8*L)
#mu = L*np.sqrt(k**2 + l**2) #dimensionless parameter governing vertical wave struture

#See APEconversion.py scripts for maximising values of mu
mu = 1.61 #Force mu to be this value. mu_crit = 2.399, mu_max (for narrow channel, l=pi/2L) = 2.01, mu_max (for wide channel, l=0) = 1.61
k = np.sqrt(((mu**2)/(L*L)) - l*l) #Choose k based on chosen value of mu

#Values from STDOUT file for MITgcm channel simulation
alpha = 2e-4 #(K^-1) thermal expansion coefficient
g = 9.81 #(ms^(-2)) gravitational acceleration
rho_0 = 999.8 #(kg/m^3) density of seawater
c_p = 3994 #(J/kg/K) specific heat capacity of seawater

#Phase speed (ms^(-1)) (taking positive root for growth), Vallis 6.86
c = U_0/2. + (U_0/mu)*cmath.sqrt((mu/2. - (1./np.tanh(mu/2.)))*(mu/2. - np.tanh(mu/2.)))

sigma = k*c.imag #(s^-1) growth rate
#sigma = 0

#Information to plot background theta distribution
theta0 = 278. #Reference temperature value, about 5 degrees celsius
dthetadz = (N**2)/(g*alpha) #Gill 6.17.24, buoyancy frequency definition
dthetady = (-f0*Lambda)/(g*alpha) #Gill 7.7.10, thermal wind


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

#Background theta distribution
def theta(y,z):
    return theta0 + dthetady*L*y + dthetadz*H*z

#Background half-slope theta distribution
def halftheta(y,z):
    return dthetady*L*y + 2*dthetadz*H*z

#Theta perturbation, Gill 13.2.5
def thetaprime(x,y,z,t):
    return (f0/(alpha*g))*(dphi_rdz(z)*np.cos(k*(L*x - c.real*T*t)) - dphi_idz(z)*np.sin(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t)


#Note that as vprime is used in the wprime expression, vprime should be decreased accordingly so that
#the amplitude of wprime is scaled correctly

#Define nondimensional velocity perturbations (see notebook for calculations)
def uprime(x,y,z,t):
    #return 100*(T/L)*l*(phi_r(z)*np.cos(k*(L*x - c.real*T*t)) - phi_i(z)*np.sin(k*(L*x - c.real*T*t)))*np.sin(l*L*y)*np.exp(sigma*T*t)
    return (T/L)*l*(phi_r(z)*np.cos(k*(L*x - c.real*T*t)) - phi_i(z)*np.sin(k*(L*x - c.real*T*t)))*np.sin(l*L*y)*np.exp(sigma*T*t)

def vprime(x,y,z,t):
    #return -100*(T/L)*k*(phi_r(z)*np.sin(k*(L*x - c.real*T*t)) + phi_i(z)*np.cos(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t)
    return -(T/L)*k*(phi_r(z)*np.sin(k*(L*x - c.real*T*t)) + phi_i(z)*np.cos(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t)

def wprime(x,y,z,t):
    #return -100*(T/H)*(f0/(N**2))*(((k*(c.real - Lambda*H*z)*dphi_rdz(z) - sigma*dphi_idz(z))*np.sin(k*(L*x - c.real*T*t)) + (k*(c.real - Lambda*H*z)*dphi_idz(z) + sigma*dphi_rdz(z))*np.cos(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t) - Lambda*(L/T)*0.01*vprime(x,y,z,t))
    return -(T/H)*(f0/(N**2))*(((k*(c.real - Lambda*H*z)*dphi_rdz(z) - sigma*dphi_idz(z))*np.sin(k*(L*x - c.real*T*t)) + (k*(c.real - Lambda*H*z)*dphi_idz(z) + sigma*dphi_rdz(z))*np.cos(k*(L*x - c.real*T*t)))*np.cos(l*L*y)*np.exp(sigma*T*t) - Lambda*(L/T)*vprime(x,y,z,t))


def WForiginvprime(y,z,t):
    return -(T/L)*k*phi_i(z)*np.cos(l*L*y)*np.exp(sigma*T*t)

def WForiginwprime(y,z,t):
    return -(T/H)*(f0/(N**2))*((k*(c.real - Lambda*H*z)*dphi_idz(z) + sigma*dphi_rdz(z))*np.cos(l*L*y)*np.exp(sigma*T*t) - Lambda*(L/T)*WForiginvprime(y,z,t))

###############################################################################################

print 'U_0 is:', U_0
print 'Lambda is:', Lambda

print '\nH is', H
print 'L (deformation radius) (m) is', L
print 'T is (seconds)', T
print 'T is (days)', T/(86400.)

print '\nk (/m) is', k
print 'l (/m) is', l
print 'mu is', mu

print '\nc (m/s) is:', c
print 'c^2 is:', np.abs(c)**2
print 'Real part of c is:', c.real
print 'sigma is:', sigma

print '\ninitial uprime amplitude (m/s) is:', l
print 'initial vprime (m/s) amplitude is:', k
print 'initial wprime (m/s) amplitude is:', -(f0/(N**2))*Lambda*k

print '\ntime (s) when vprime amplitude is equal to U_0 is:', (1./sigma)*np.log(U_0/(k))
print 'in periods of T, that is:', ((1./sigma)*np.log(U_0/(k)))/T

print '\nMax Eady growth rate (/s) is approximately:', (-0.31*Lambda*f0)/(N)
print 'Which corresponds to a timescale (s) of:', N/(-0.31*Lambda*f0)
print 'in days, that is:', N/(-0.31*Lambda*f0*86400.)
print 'in periods of T, that is:', N/(-0.31*Lambda*f0*T)

print '\nLengthscale of maximum instability (3.9*L) (m) is approximately:', 3.9*L

print '\nRossby number is:', U_0/np.abs(f0*L)
print 'Richardson number is:', (N**2)/(Lambda**2)
