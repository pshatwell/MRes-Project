import numpy as np

################################# LIST OF PHYSICAL PARAMETERS #############################################################

f0 = 1e-4 #(s^-1) f-plane approx. for Coriolis parameter
N = 1e-2 #(s^-1) buoyancy frequency
k = 1e-6 #(m^-1) zonal wavenumber of growing wave
H_R = np.divide(f0, N*k) #(m) Rossby height
H = 1e4 #(m) positions of boundaries in Eady problem at H and -H
Hratio = np.divide(H, H_R)
deltaU = 10 #(ms^(-1)) velocity difference between boundaries at H and -H
shear = np.divide(deltaU, 2*H) #(s^-1) velocity shear
sigma_max = 0.3098*np.divide(f0,N)*shear #maximum growth rate
