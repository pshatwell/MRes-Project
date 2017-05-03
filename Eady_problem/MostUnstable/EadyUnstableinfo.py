import numpy as np

#Parameters and functions for EadyUnstable scripts

################################# LIST OF PHYSICAL PARAMETERS #############################################################

#PARAMETERS ARE DIFFERENT NOW (6/4/17). TAKEN FROM VALLIS BOOK (P. 268) FOR MOST UNSTABLE CASE FOR OCEAN.

#Typical parameters for ocean, for most unstable (fastest growing) Eady wave

f0 = -1e-4 #(s^-1) f-plane approx. for Coriolis parameter (southern hemisphere)
#N = 1e-3 #(s^-1) buoyancy frequency OLD
N = 1e-2 #(s^-1) buoyancy frequency
#L = 1e4 #(m) typical length scale OLD
L = 1e5 #(m) typical length scale of deformation radius in ocean
#k = 1e-4 #(m^-1) zonal wavenumber of growing wave OLD
k = 0.8031/L #(m^-1) zonal wavenumber for most unstable growing wave, Gill 13.3.12
l = np.pi/(2.*L) #(m^-1) meridional wavenumber (= 0 for most unstable case), defining zonal channel

H_R = np.abs(f0/(N*k)) #(m) Rossby height
H = 0.8031*H_R #(m) positions of boundaries at H and -H for most unstable wave, Gill 13.3.12
Hratio = np.abs(H/H_R)
#U = 0.25 #(ms^(-1)) mean flow zonal velocity magnitude at boundaries H and -H OLD
U = 0.1 #(ms^(-1)) typical mean flow zonal velocity magnitude at boundaries H and -H
shear = U/H #(s^-1) velocity shear
#T = 10.*np.abs(1./f0) #(s) typical time scale ~ a day OLD
T = L/U #Eady timescale
sigma_max = np.abs(0.3098*(f0/N)*shear) #maximum growth rate
#sigma_max = 0
alpha = 2e-4 #(K^-1) thermal expansion coefficient
g = 9.8 #(ms^(-2)) gravitational acceleration

c = U/L #phase speed, dimensions of s^(-1)

#Information to attempt to plot background theta distribution
theta0 = 280. #Not important at the moment
dthetadz = (theta0*(N**2))/g
dthetady = -(f0*shear)/(g*alpha)

################################# FUNCTIONS #############################################################

#Define streamfunction and its derivatives as velocity perturbations
#Make x, z and t all non-dimensional

#CHECK TIME DERIVATIVES FOR TRANSFORMATION, TRY CHANGING AMPLITUDES

#ADDED COSINE FUNCTION IN y TO PHIPRIME. NEW EXPRESSIONS BENEATH OLD. OLD EXPRESSIONS HAVE NO y-DEPENDENCE.

def phiprime(x,y,z,t): #Streamfunction perturbation for fastest growing Eady mode, Gill 13.3.15
    #return ((np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))
    return ((np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))*np.cos(l*L*y)

def phiprime_x(x,y,z,t):
    #return (k*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))
    return (k*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))*np.cos(l*L*y)

def phiprime_y(x,y,z,t):
    return (-l*(np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t)*np.sin(l*L*y))

def phiprime_z(x,y,z,t):
    #return ((1./H_R)*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))
    return ((1./H_R)*(np.cos(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))) + np.sin(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))))*np.exp(sigma_max*T*t))*np.cos(l*L*y)

def phiprime_t(x,y,z,t):
    #return (sigma_max*phiprime(x,z,t) - c*phiprime_x(x,z,t))
    return sigma_max*phiprime(x,y,z,t)

def phiprime_zx(x,y,z,t):
    #return ((k/H_R)*(np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))
    return ((k/H_R)*(np.cos(k*L*(x-c*T*t))*((np.sinh(z*Hratio))/(np.cosh(Hratio))) - np.sin(k*L*(x-c*T*t))*((np.cosh(z*Hratio))/(np.sinh(Hratio))))*np.exp(sigma_max*T*t))*np.cos(l*L*y)

def phiprime_zt(x,y,z,t):
    #return (sigma_max*phiprime_z(x,z,t) - c*phiprime_zx(x,z,t))
    return sigma_max*phiprime_z(x,y,z,t)

def umeanflow(z): #Non-dimensional zonal velocity mean shear flow
    return (1./(L*np.abs(f0)))*shear*H*z

def uprime(x,y,z,t):
    return (1./(L*np.abs(f0)))*phiprime_y(x,y,z,t)

def vprime(x,y,z,t): #Non-dimensional meridional velocity perturbation, Gill 12.9.3
    return (1./(L*(f0*f0)))*phiprime_x(x,y,z,t)

#CHECK THE SIGN OF f0 HERE??
#AT THE MOMENT, ISOTHERMALS SLOPE IN THE WRONG DIRECTION??
def wprime(x,y,z,t): #Non-dimensional vertical velocity perturbation, Gill 12.9.6
    return (1./(H*np.abs(f0)))*(1./N**2)*(shear*phiprime_x(x,y,z,t) - phiprime_zt(x,y,z,t) - shear*H*z*phiprime_zx(x,y,z,t))

def thetaprime(x,y,z,t): #Potential temperature perturbation, Gill 13.2.5
    return (1./(alpha*g))*phiprime_z(x,y,z,t)

def theta(y,z): #Background theta distribution
    return theta0 + dthetady*y + dthetadz*z
