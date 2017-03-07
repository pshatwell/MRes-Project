import numpy as np

#Parameters and functions for EadyUnstable script. In the Earth frame.
#AT THE MOMENT, THINGS PROPAGATE WESTWARDS. IS THIS CORRECT?

#CAN YOU JUST CHANGE VARIABLES TO x-cTt? WHAT ABOUT NEW TIME DERIVATIVES?

################################# LIST OF PHYSICAL PARAMETERS #############################################################

#Typical parameters for ocean, for most unstable (fastest growing) Eady wave

L = 1e4 #(m) typical length scale
f0 = -1e-4 #(s^-1) f-plane approx. for Coriolis parameter (southern hemisphere)
T = 10*np.abs(1./f0) #(s) typical time scale ~ a day
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
c = 0.0001 #(ms^(-1)) phase speed of wave CHECK WHAT A PHYSICAL VALUE SHOULD BE

################################# FUNCTIONS #############################################################

#Define streamfunction and its derivatives as velocity perturbations
#Make x, z and t all non-dimensional

def phiprime(x,z,t): #Streamfunction perturbation for fastest growing Eady mode, Gill 13.3.15
    return ((np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))

def phiprime_x(x,z,t):
    return (k*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))

def phiprime_z(x,z,t):
    return ((1./H_R)*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))

def phiprime_t(x,z,t):
    #return (sigma_max*phiprime(x,z,t) - c*phiprime_x(x,z,t))
    return sigma_max*phiprime(x,z,t)

def phiprime_zx(x,z,t):
    return ((k/H_R)*(np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))

def phiprime_zt(x,z,t):
    #return (sigma_max*phiprime_z(x,z,t) - c*phiprime_zx(x,z,t))
    return sigma_max*phiprime_z(x,z,t)

def vprime(x,z,t): #Non-dimensional meridional velocity perturbation, Gill 12.9.3
    return (1./(L*(f0**2)))*phiprime_x(x,z,t)

def wprime(x,z,t): #Non-dimensional vertical velocity perturbation, Gill 12.9.6
    return (1./(H*f0))*(1./N**2)*(shear*phiprime_x(x,z,t) - phiprime_zt(x,z,t) - shear*H*z*phiprime_zx(x,z,t))

def thetaprime(x,z,t): #Potential temperature perturbation, Gill 13.2.5
    return (1./(alpha*g))*phiprime_z(x,z,t)
