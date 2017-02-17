import numpy as np
from sympy import *
from RossbyInfo import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

###############################################################################################

#PRINT OUT PHYSICAL QUANTITIES FOR THE WAVE AND DIMENSIONAL SCALINGS

print 'phase speed, c (ms^(-1)) is:', c
print 'horizontal wavelength (m) is:', wavelength

print 'amplitude of psiprime (in m^2*s^(-1)) is:', psi0
print 'amplitude of uprime (in ms^(-1)), psi0*pi/L, is:', (psi0*np.pi)/L #amplitude of uprime
print 'amplitude of vprime (in ms^(-1)), psi0*k, is:', psi0*k #amplitude of vprime

print 'k is:', k
print 'L is:', L
print 'k*L is:', k*L

print 'c is:', c
print 'T is:', T
print 'k*c*T is:', k*c*T

###############################################################################################

def main():
    #Define dimensions of arrays
    xvalues = np.linspace(0,2,200)
    yvalues = np.linspace(0,1,100) #Channel width in y direction is 1*L

    time = 0 #Choose time step (in units of T)

    xlength = len(xvalues)
    ylength = len(yvalues)

###############################################################################################

    #CREATE ARRAY FOR PSIPRIME VALUES
    psiprime_matrix = np.zeros((ylength,xlength)) #Note x and y have switched so plotted output is oriented correctly

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            psiprime_matrix[i,j] = psiprime(x=xvalues[j], y=yvalues[i], t=time)

    print 'type of psiprime_matrix is:', type(psiprime_matrix)
    print 'shape of psiprime_matrix is:', psiprime_matrix.shape

###############################################################################################

    #CREATE ARRAY FOR UPRIME VALUES
    uprime_matrix = np.zeros((ylength,xlength)) #Note x and y have switched so plotted output is oriented correctly

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            uprime_matrix[i,j] = uprime(x=xvalues[j], y=yvalues[i], t=time)

    print 'type of uprime_matrix is:', type(uprime_matrix)
    print 'shape of uprime_matrix is:', uprime_matrix.shape

###############################################################################################

    #CREATE ARRAY FOR VPRIME VALUES
    vprime_matrix = np.zeros((ylength,xlength))

    for i in range(0, ylength, 1):
        for j in range(0, xlength, 1):
            vprime_matrix[i,j] = vprime(x=xvalues[j], y=yvalues[i], t=time)


    print 'type of vprime_matrix is:', type(vprime_matrix)
    print 'shape of vprime_matrix is:', vprime_matrix.shape

    #Check boundary conditions for v
    print 'vprime at y = 0 is:', vprime_matrix[0,10]
    print 'vprime at y = L is:', vprime_matrix[-1,10]

###############################################################################################

    #DEFINE ODE TO INTEGRATE FOR PARCEL TRAJECTORIES
    def velocity(s,t):
        x, y = s
        dsdt = [uprime(x=x,y=y,t=t), vprime(x=x,y=y,t=t)]
        #dsdt = [testuprime(x=x,y=y,t=t), testvprime(x=x,y=y,t=t)]
        return dsdt

    t=np.linspace(0,100,200)

    #s0_a = [0.47, 0.1] #oscillatory trajectory
    #s0_e = [0.7857, 0.5] #oscillatory trajectory

    p=0.2 #Choose point to determine cluster of initial points for trajectories

    #Create cluster of 9 points to evolve in time
    s0_a = [p, p]
    s0_b = [p, p+0.1]
    s0_c = [p, p+0.2]
    s0_d = [p+0.1, p]
    s0_e = [p+0.1, p+0.1]
    s0_f = [p+0.1, p+0.2]
    s0_g = [p+0.2, p]
    s0_h = [p+0.2, p+0.1]
    s0_i = [p+0.2, p+0.2]

    sol_a = odeint(velocity, s0_a, t)
    sol_b = odeint(velocity, s0_b, t)
    sol_c = odeint(velocity, s0_c, t)
    sol_d = odeint(velocity, s0_d, t)
    sol_e = odeint(velocity, s0_e, t)
    sol_f = odeint(velocity, s0_f, t)
    sol_g = odeint(velocity, s0_g, t)
    sol_h = odeint(velocity, s0_h, t)
    sol_i = odeint(velocity, s0_i, t)

    print 'type of trajectory solution is:', type(sol_a)
    print 'shape of trajectory solution is:', sol_a.shape


###############################################################################################

    #TRANSFORM TRAJECTORIES TO WAVE REFERENCE FRAME
    shift = np.zeros_like(sol_a)
    shift[:,0] = c*t

    print 'shape of sol_a is:', sol_a.shape
    print 'shape of shift is:', shift.shape

    rel_sol_a = sol_a - shift #new parcel trajectory relative to Rossby wave

###############################################################################################

    #DO ALL THE PLOTTING

    fig1 = plt.figure()
    plt.set_cmap('inferno') #Choose colour map

    #Note 'extents' match the sizes of xvalues and yvalues defined at the beginning.
    #(In units of L)
    ax1=fig1.add_subplot(311)
    #ax1.set_title('Streamfunction')
    psi_show = ax1.imshow(psiprime_matrix[:,:], origin='lower', extent=[0,2,0,1], aspect='auto')
    plt.colorbar(psi_show) #Note units are m^2*s^(-1)
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')

    ax2 = fig1.add_subplot(312)
    #ax2.set_title('u')
    u_show = ax2.imshow(uprime_matrix[:,:], origin='lower', extent=[0,2,0,1], aspect='auto')
    plt.colorbar(u_show) #Note units are m*s^(-1)
    ax2.set_xlabel('x (L)')
    ax2.set_ylabel('y (L)')

    ax3 = fig1.add_subplot(313)
    #ax3.set_title('v')
    v_show = ax3.imshow(vprime_matrix[:,:], origin='lower', extent=[0,2,0,1], aspect='auto')
    plt.colorbar(v_show) #Note units are m*s^(-1)
    ax3.set_xlabel('x (L)')
    ax3.set_ylabel('y (L)')


    fig2 = plt.figure()
    ax4 = fig2.add_subplot(111)
    ax4.set_xlabel('x (L)')
    ax4.set_ylabel('y (L)')

    ax4.plot(sol_a[:,0], sol_a[:,1], linewidth=1.5, label='a')
    ax4.plot(sol_b[:,0], sol_b[:,1], linewidth=1.5, label='b')
    ax4.plot(sol_c[:,0], sol_c[:,1], linewidth=1.5, label='c')
    ax4.plot(sol_d[:,0], sol_d[:,1], linewidth=1.5, label='d')
    ax4.plot(sol_e[:,0], sol_e[:,1], linewidth=1.5, label='e')
    ax4.plot(sol_f[:,0], sol_f[:,1], linewidth=1.5, label='f')
    ax4.plot(sol_g[:,0], sol_g[:,1], linewidth=1.5, label='g')
    ax4.plot(sol_h[:,0], sol_h[:,1], linewidth=1.5, label='h')
    ax4.plot(sol_i[:,0], sol_i[:,1], linewidth=1.5, label='i')

    ax4.plot(rel_sol_a[:,0], rel_sol_a[:,1], color='black', label ='a in RW frame')

    #Overlay initial streamfunction to see if trajectories make sense.
    ax4.imshow(psiprime_matrix[:,:], origin='lower', extent=[0,2,0,1], aspect='auto')

    #plt.legend(loc='upper left', bbox_to_anchor=(0, 1))

    plt.show()

main()
