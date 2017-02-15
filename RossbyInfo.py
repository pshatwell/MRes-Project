import numpy as np

#PARAMETER VALUES AND FUNCTIONS USED FOR ROSSBY WAVE SCRIPTS
#SET DIMENSIONAL SCALES FOR VARIABLES

L = 1e5 #(m) typical length scale
f0 = 1e-4 #(s^-1) Coriolis parameter
T = 1./f0 #(s) typical time scale
#T = 1000
k = 2e-5 #(m^-1) horizontal wavenumber
a = 3e4 #(m) deformation radius
psi0 = 3e3 #(m^2*s^-1) wave amplitude
beta = 3e-11 #(m^-1*s^-1) latitudinal variation of Coriolis parameter
c = -beta*(a**2) #(m*s^-1) phase speed of wave

A_und = (psi0*np.pi)/((L**2)*f0) #non-dimensional amplitude of uprime
A_vnd = (psi0*k)/(L*f0) #non-dimensional amplitude of vprime

#Define streamfunction and its derivatives as velocity perturbations
#Make x,y and t all non-dimensional
def psiprime(x,y,t):
    return psi0*np.cos(k*L*x - c*T*t)*np.sin(np.pi*y)

def uprime(x,y,t): #zonal velocity, equal to -psiprime_y
    #return -A_und*np.cos(k*L*x - c*T*t)*np.cos(np.pi*y)
    return -np.cos(k*L*x - c*T*t)*np.cos(np.pi*y)

def vprime(x,y,t): #meridional velocity, equal to psiprime_x
    #return -A_vnd*np.sin(k*L*x - c*T*t)*np.sin(np.pi*y)
    return -np.sin(k*L*x - c*T*t)*np.sin(np.pi*y)

def testuprime(x,y,t):
    return -np.cos(2*x - t)*np.cos(2*y)

def testvprime(x,y,t):
    return -np.sin(2*x - t)*np.sin(2*y)
