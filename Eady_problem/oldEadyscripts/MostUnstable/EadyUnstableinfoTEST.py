import numpy as np
import cmath

#Parameters and functions for EadyUnstable scripts

################################# LIST OF PHYSICAL PARAMETERS #############################################################

#Typical parameters for ocean, for most unstable (fastest growing) Eady wave (from Vallis p. 268)

f0 = -1e-4 #(s^-1) f-plane approx. for Coriolis parameter (southern hemisphere)
N = 2e-3 #(s^-1) buoyancy frequency (for SO)
L = 1e5 #(m) typical length scale of deformation radius in ocean
k = 1e-5 #(m^-1) zonal wavenumber (=0.8031/L for most unstable case, Gill 13.3.12)
l = np.pi/(2.*L) #(m^-1) meridional wavenumber (= 0 for most unstable case), defining zonal channel
#l=0
kappa = np.sqrt(k**2 + l**2) #(m^-1) total wavenumber

H_R = np.abs(f0/(N*kappa)) #(m) Rossby height
H = 15e2 #(m) height of boundaries (to fit 3km depth from MITgcm simulation comparison)
Hratio = H/H_R
U0 = 0.1 #(ms^(-1)) typical mean flow zonal velocity magnitude at boundaries H and -H
shear = U0/H #(s^-1) velocity shear
T = L/U0 #(s) Eady timescale
alpha = 2e-4 #(K^-1) thermal expansion coefficient
g = 9.8 #(ms^(-2)) gravitational acceleration

#Phase speed (in units of s^-1), Gill 13.3.5 + U0/L
c = (U0 + shear*cmath.sqrt((H*np.tanh(Hratio) - H_R)*(H*(1./np.tanh(Hratio)) - H_R)))/L

sigma = k*c.imag*L #(s^-1) growth rate

#Information to plot background theta distribution
theta0 = 280. #Reference temperature value
dthetadz = (N**2)/(g*alpha) #Gill 6.17.24, buoyancy frequency definition
dthetady = (-f0*shear)/(g*alpha) #Gill 7.7.10, thermal wind

################################# FUNCTIONS #############################################################

#Define streamfunction and its derivatives as velocity perturbations
#Make x, z and t all non-dimensional

def phiprime(x,y,z,t): #Streamfunction (geopotential) perturbation, Gill 13.3
    #return ((np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))
    #return ((np.cos(k*L*(x-c.real*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c.real*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma*T*t))*np.cos(l*L*y)
    return (np.cos(k*L*(x - c.real*T*t))*np.sinh(z*Hratio) + np.sin(k*L*(x - c.real*T*t))*np.cosh(z*Hratio))*np.cos(l*L*y)*np.exp(sigma*T*t)

def phiprime_x(x,y,z,t):
    #return (k*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))
    #return (k*(np.cos(k*L*(x-c.real*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c.real*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma*T*t))*np.cos(l*L*y)
    return -k*(np.sin(k*L*(x - c.real*T*t))*np.sinh(z*Hratio) - np.cos(k*L*(x - c.real*T*t))*np.cosh(z*Hratio))*np.cos(l*L*y)*np.exp(sigma*T*t)

def phiprime_y(x,y,z,t):
    #return (-l*(np.cos(k*L*(x-c.real*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c.real*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma*T*t)*np.sin(l*L*y))
    return -l*(np.cos(k*L*(x - c.real*T*t))*np.sinh(z*Hratio) + np.sin(k*L*(x - c.real*T*t))*np.cosh(z*Hratio))*np.sin(l*L*y)*np.exp(sigma*T*t)

def phiprime_z(x,y,z,t):
    #return ((1./H_R)*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))
    #return ((1./H_R)*(np.cos(k*L*(x-c.real*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c.real*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma*T*t))*np.cos(l*L*y)
    return (1./H_R)*(np.cos(k*L*(x - c.real*T*t))*np.cosh(z*Hratio) + np.sin(k*L*(x - c.real*T*t))*np.sinh(z*Hratio))*np.cos(l*L*y)*np.exp(sigma*T*t)

def phiprime_t(x,y,z,t):
    #return (sigma_max*phiprime(x,z,t) - c*phiprime_x(x,z,t))
    return sigma*phiprime(x,y,z,t)

def phiprime_zx(x,y,z,t):
    #return ((k/H_R)*(np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))
    #return ((k/H_R)*(np.cos(k*L*(x-c.real*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c.real*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma*T*t))*np.cos(l*L*y)
    return (-k/H_R)*(np.sin(k*L*(x - c.real*T*t))*np.cosh(z*Hratio) - np.cos(k*L*(x - c.real*T*t))*np.sinh(z*Hratio))*np.cos(l*L*y)*np.exp(sigma*T*t)

def phiprime_zt(x,y,z,t):
    #return (sigma_max*phiprime_z(x,z,t) - c*phiprime_zx(x,z,t))
    return sigma*phiprime_z(x,y,z,t)

def umeanflow(z): #Non-dimensional zonal velocity mean shear flow
    return (T/L)*shear*H*z

def uprime(x,y,z,t): #Non-dimensional zonal velocity perturbation
    return (T/L)*(-1./f0)*phiprime_y(x,y,z,t)

def vprime(x,y,z,t): #Non-dimensional meridional velocity perturbation, Gill 12.9.3
    return (T/L)*(1./f0)*phiprime_x(x,y,z,t)

def wprime(x,y,z,t): #Non-dimensional vertical velocity perturbation, Gill 12.9.6
    return (T/H)*(1./N**2)*(shear*phiprime_x(x,y,z,t) - phiprime_zt(x,y,z,t) - shear*H*z*phiprime_zx(x,y,z,t))

def thetaprime(x,y,z,t): #Potential temperature perturbation, Gill 13.2.5
    return (1./(alpha*g))*phiprime_z(x,y,z,t)

def theta(y,z): #Background theta distribution
    return theta0 + dthetady*y + dthetadz*z
