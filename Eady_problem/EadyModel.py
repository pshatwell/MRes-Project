import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from Eadyinfo import *
from mpl_toolkits.mplot3d import Axes3D

def main():
    xvalues = np.linspace(0,6,50)
    yvalues = np.linspace(0,6,50)
    zvalues = np.linspace(-1,1,50)

    time = 10

    xlength = len(xvalues)
    ylength = len(yvalues)
    zlength = len(zvalues)

    #Create empty matrix for streamfunction
    phi_matrix = np.zeros((zlength,xlength,ylength))

    #Evaluate streamfunction matrix values
    for i in range(0, zlength, 1):
        for j in range(0, xlength, 1):
            for m in range(0, ylength, 1):
                phi_matrix[i,j,m] = streamfunction(x=xvalues[j], y = yvalues[m], z=zvalues[i], t=time)

    print 'type of phi_matrix is:', type(phi_matrix)
    print 'shape of phi_matrix is:', phi_matrix.shape
###############################################################################################

    fig = plt.figure()
    plt.set_cmap('inferno')
    ax1 = fig.add_subplot(211)
    ax1.set_title('Streamfunction')
    ax1.set_xlabel('x (L)')
    ax1.set_ylabel('y (L)')
    ax1.contourf(phi_matrix[0,:,:], origin='lower', aspect='auto')

    ax2 = fig.add_subplot(212)
    ax2.set_title('Streamfunction')
    ax2.set_xlabel('y (L)')
    ax2.set_ylabel('z (L)')
    ax2.contourf(phi_matrix[:,0,:], origin='lower', aspect='auto')

    plt.show()

main()
