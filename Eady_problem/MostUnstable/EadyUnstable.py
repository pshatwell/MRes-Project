import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint


#Choose file to import to give trajectories in Earth frame or Wave frame

#from EadyUnstableinfo_Waveframe import *
from EadyUnstableinfo_Earthframe import *

###############################################################################################

print 'L is (m):', L
print 'H is (m):', H
print 'H_R is (m):', H_R
print 'H/H_R is:', Hratio
print 'T is (s):', T
print 'k is (m^-1):', k
print 'sigma_max is (s^-1):', sigma_max
print 'e-folding time (T) is (sigma_max*T)^(-1):', 1./(sigma_max*T)
print 'N is (s^-1):', N
print 'U is (m*s^-1):', U
print 'velocity shear is (s^-1):', shear
print 'velocity scaling (to dimensionalise) is:', L*f0
print 'phase speed c is (m*s^-1):', c
print 'angular frequency, k*c is (s^-1):', k*c

###############################################################################################

#xslice defines meridional plane to view trajectories in
#parameter p defines position of initial cluster

def main(p, xslice):
    xvalues = np.linspace(0,6,100)
    zvalues = np.linspace(-1,1,100)

    time = 0

    xlength = len(xvalues)
    zlength = len(zvalues)

###############################################################################################

    #Create empty matrix for streamfunction
    phiprime_matrix = np.zeros((zlength,xlength))

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            phiprime_matrix[i,j] = phiprime(x=xvalues[j], z=zvalues[i], t=time)

    #print 'type of phiprime_matrix is:', type(phiprime_matrix)
    #print 'shape of phiprime_matrix is:', phiprime_matrix.shape


    #Create empty matrix for meridional velocity perturbation
    vprime_matrix = np.zeros((zlength,xlength))

    #Evaluate vprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            vprime_matrix[i,j] = vprime(x=xvalues[j], z=zvalues[i], t=time)

    #print 'type of vprime_matrix is:', type(vprime_matrix)
    #print 'shape of vprime_matrix is:', vprime_matrix.shape


    #Create empty matrix for vertical velocity perturbation
    wprime_matrix = np.zeros((zlength,xlength))

    #Evaluate wprime matrix values
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            wprime_matrix[i,j] = wprime(x=xvalues[j], z=zvalues[i], t=time)

    #print 'type of wprime_matrix is:', type(wprime_matrix)
    #print 'shape of wprime_matrix is:', wprime_matrix.shape

    print 'wprime at z=-H is:', wprime_matrix[0,10]
    print 'wprime at z=H is:', wprime_matrix[-1,10]

    #Create empty matrix for potential temperature perturbation
    thetaprime_matrix = np.zeros((zlength,xlength))

    #Evaluate thetaprime matrix values:
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            thetaprime_matrix[i,j] = thetaprime(x=xvalues[j], z=zvalues[i], t=time)


###############################################################################################

    def velocity(s,t):
        y,z = s
        dsdt = [vprime(x=xslice,z=z,t=t), wprime(x=xslice,z=z,t=t)]
        #dsdt = [vprime(x=4,z=z,t=t), 0.001]
        return dsdt

    t = np.linspace(0,10,200)

    s0_a = [p, p/5.]
    #s0_b = [p-0.1, (p-0.1)/5.]
    #s0_c = [p+0.1, (p+0.1)/5.]
    #s0_d = [p-0.2, (p-0.2)/5.]
    #s0_e = [p+0.2, (p+0.2)/5.]
    #s0_f = [p-0.3, (p-0.3)/5.]
    #s0_g = [p+0.3, (p+0.3)/5.]
    #s0_h = [p-0.4, (p-0.4)/5.]
    #s0_i = [p+0.4, (p+0.4)/5.]

    sol_a = odeint(velocity, s0_a, t)
    #sol_b = odeint(velocity, s0_b, t)
    #sol_c = odeint(velocity, s0_c, t)
    #sol_d = odeint(velocity, s0_d, t)
    #sol_e = odeint(velocity, s0_e, t)
    #sol_f = odeint(velocity, s0_f, t)
    #sol_g = odeint(velocity, s0_g, t)
    #sol_h = odeint(velocity, s0_h, t)
    #sol_i = odeint(velocity, s0_i, t)

###############################################################################################

    fig = plt.figure()
    plt.set_cmap('inferno')
    ax1 = fig.add_subplot(221)
    ax1.set_title('Streamfunction')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('z (H)')
    phicontour = ax1.contourf(phiprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(phicontour)

    ax2 = fig.add_subplot(222)
    ax2.set_title('v')
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('z (H)')
    vcontour = ax2.contourf(vprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(vcontour)

    ax3 = fig.add_subplot(223)
    ax3.set_title('w')
    ax3.set_xlabel('x (L)')
    ax3.set_ylabel('z (H)')
    wcontour = ax3.contourf(wprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(wcontour)

    ax4 = fig.add_subplot(224)
    ax4.set_title('theta')
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('z (H)')
    thetacontour = ax4.contourf(thetaprime_matrix, origin='lower', extent=[0,6,-1,1], aspect='auto')
    plt.colorbar(thetacontour)

    #plt.savefig('figures/EadyUnstableSolution.pdf')

    fig2 = plt.figure()
    ax5 = fig2.add_subplot(111)
    ax5.set_xlabel('y (L)')
    ax5.set_ylabel('z (H)')
    ax5.plot(sol_a[:,0],sol_a[:,1])
    #ax5.plot(sol_b[:,0],sol_b[:,1])
    #ax5.plot(sol_c[:,0],sol_c[:,1])
    #ax5.plot(sol_d[:,0],sol_d[:,1])
    #ax5.plot(sol_e[:,0],sol_e[:,1])
    #ax5.plot(sol_f[:,0],sol_f[:,1])
    #ax5.plot(sol_g[:,0],sol_g[:,1])
    #ax5.plot(sol_h[:,0],sol_h[:,1])
    #ax5.plot(sol_i[:,0],sol_i[:,1])

    plt.savefig('figures/c_%s_p_%s_xslice_%s.pdf' % (str(c), str(p), str(xslice)))


    plt.show()

#xslices = [i for i in np.arange(0,6,0.1)]

#for i in xslices:
#    main(p=0, xslice=i)

main(p=0, xslice=0)
