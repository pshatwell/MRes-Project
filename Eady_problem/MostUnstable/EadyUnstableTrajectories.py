import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

from EadyUnstableinfo import *
from EadyUnstableSolution import *

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

#Use main function to avoid global variables
def main(p, xslice):

###############################################################################################

    #Define velocity functions for parcel trajectories
    def velocity(s,t):
        y,z = s
        dsdt = [vprime(x=xslice,z=z,t=t), wprime(x=xslice,z=z,t=t)]
        return dsdt

    def velocity3d(s,t):
        x,y,z = s
        dsdt = [umeanflow(z=z), vprime(x=x,z=z,t=t), wprime(x=x,z=z,t=t)]
        return dsdt

###############################################################################################

    tmin = 0
    tmax = 5

    t = np.linspace(tmin, tmax, 200)

    #THESE TRAJECTORIES MAY JUST BE NONSENSE. ISN'T THE MEAN FLOW IMPORTANT??
    #Solve for parcel trajectories
    s0_a = [p, p/5.]
    s0_b = [p-0.1, (p-0.1)/5.]
    s0_c = [p+0.1, (p+0.1)/5.]
    s0_d = [p-0.2, (p-0.2)/5.]
    s0_e = [p+0.2, (p+0.2)/5.]
    s0_f = [p-0.3, (p-0.3)/5.]
    s0_g = [p+0.3, (p+0.3)/5.]
    s0_h = [p-0.4, (p-0.4)/5.]
    s0_i = [p+0.4, (p+0.4)/5.]

    sol_a = odeint(velocity, s0_a, t)
    sol_b = odeint(velocity, s0_b, t)
    sol_c = odeint(velocity, s0_c, t)
    sol_d = odeint(velocity, s0_d, t)
    sol_e = odeint(velocity, s0_e, t)
    sol_f = odeint(velocity, s0_f, t)
    sol_g = odeint(velocity, s0_g, t)
    sol_h = odeint(velocity, s0_h, t)
    sol_i = odeint(velocity, s0_i, t)

###############################################################################################

    #Solve for parcel trajectories in 3d
    s0_a_3d = [p, p, p/5.]
    s0_b_3d = [p, p-0.1, (p-0.1)/5.]
    s0_c_3d = [p, p+0.1, (p+0.1)/5.]
    s0_d_3d = [p, p-0.2, (p-0.2)/5.]
    s0_e_3d = [p, p+0.2, (p+0.2)/5.]
    s0_f_3d = [p, p-0.3, (p-0.3)/5.]
    s0_g_3d = [p, p+0.3, (p+0.3)/5.]
    s0_h_3d = [p, p-0.4, (p-0.4)/5.]
    s0_i_3d = [p, p+0.4, (p+0.4)/5.]

    sol_a_3d = odeint(velocity3d, s0_a_3d, t)
    sol_b_3d = odeint(velocity3d, s0_b_3d, t)
    sol_c_3d = odeint(velocity3d, s0_c_3d, t)
    sol_d_3d = odeint(velocity3d, s0_d_3d, t)
    sol_e_3d = odeint(velocity3d, s0_e_3d, t)
    sol_f_3d = odeint(velocity3d, s0_f_3d, t)
    sol_g_3d = odeint(velocity3d, s0_g_3d, t)
    sol_h_3d = odeint(velocity3d, s0_h_3d, t)
    sol_i_3d = odeint(velocity3d, s0_i_3d, t)

###############################################################################################

    #Plot trajectories

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('y (L)')
    ax.set_ylabel('z (H)')
    #thetacontour = ax.contourf(theta_matrix, origin='lower', extent=[ymin, ymax, zmin, zmax], aspect='auto')
    #plt.colorbar(thetacontour)

    ax.plot(sol_a[:,0],sol_a[:,1])
    ax.plot(sol_b[:,0],sol_b[:,1])
    ax.plot(sol_c[:,0],sol_c[:,1])
    ax.plot(sol_d[:,0],sol_d[:,1])
    ax.plot(sol_e[:,0],sol_e[:,1])
    ax.plot(sol_f[:,0],sol_f[:,1])
    ax.plot(sol_g[:,0],sol_g[:,1])
    ax.plot(sol_h[:,0],sol_h[:,1])
    ax.plot(sol_i[:,0],sol_i[:,1])

    #plt.savefig('figures/EarthFrame/c_%s_p_%s_xslice_%s.pdf' % (str(c), str(p), str(xslice)))
    #plt.savefig('figures/WaveFrame/c_%s_p_%s_xslice_%s.pdf' % (str(c), str(p), str(xslice)))


    #3d plotting

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection = '3d')
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('y (L)')
    ax2.set_zlabel('z (H)')
    ax2.plot(sol_a_3d[:,0], sol_a_3d[:,1], sol_a_3d[:,2])
    ax2.plot(sol_b_3d[:,0], sol_b_3d[:,1], sol_b_3d[:,2])
    ax2.plot(sol_c_3d[:,0], sol_c_3d[:,1], sol_c_3d[:,2])
    ax2.plot(sol_d_3d[:,0], sol_d_3d[:,1], sol_d_3d[:,2])
    ax2.plot(sol_e_3d[:,0], sol_e_3d[:,1], sol_e_3d[:,2])
    ax2.plot(sol_f_3d[:,0], sol_f_3d[:,1], sol_f_3d[:,2])
    ax2.plot(sol_g_3d[:,0], sol_g_3d[:,1], sol_g_3d[:,2])
    ax2.plot(sol_h_3d[:,0], sol_h_3d[:,1], sol_h_3d[:,2])
    ax2.plot(sol_i_3d[:,0], sol_i_3d[:,1], sol_i_3d[:,2])

    plt.show()


#xslices = [i for i in np.arange(0,6,0.1)]

#for i in xslices:
#    main(p=0, xslice=i)

#pvalues = [i for i in np.arange(-4,4,0.5)]

#for i in pvalues:
    #main(p=i, xslice=0)

main(p=1.8, xslice=0)