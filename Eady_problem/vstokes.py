from __future__ import division

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

'''Plotting the variation of vertical Stokes drifts with initial z position from Wallace trajectories script'''

zposthings = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
zposthings2 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
zposthings3 = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.45,0.5,0.6,0.65,0.75,0.8,0.85,0.9,0.95]

#Data using MEAN values for vertical DRIFTS
ups3 = [0.0023,0.0064,0.0066,0.0168,0.0182,0.0339,0.0203,0.0361,0.0553,0.0422,0.0209,0.0057,0.0062,0.0037,0.0047,0.0025]
downs3 = [0.0029,0.0036,0.0045,0.0062,0.0054,0.0149,0.0333,0.0377,0.0135,0.0162,0.0456,0.0271,0.0033,0.0103,0.0051,0.0024]

#ALL DATA BELOW USING MAX AND MIN VALUES FOR VERTICAL DRIFTS
#More data for stop=31.6
ups = [0.0272,0.069,0.1357,0.1753,0.1182,0.1853,0.195,0.1499,0.1835,0.1785,0.1655,0.1482,0.1247,0.0909,0.0631,0.0385,0.0477,0.0212]
downs = [0.0234,0.0406,0.0478,0.0729,0.1042,0.1422,0.165,0.1736,0.209,0.1723,0.1671,0.2004,0.1529,0.2172,0.0976,0.1021,0.0537,0.0293]

#Data for stop=25.3
ups1 = [0.0107,0.0601,0.1406,0.0021,0.1835,0.0639,0.0142,0.0631,0.0196]
downs1 = [0.0267,0.0381,0.0374,0.1225,0.0037,0.1152,0.1528,0.0301,0.0465]

#Data for stop=31.6
ups2 = [0.069,0.1753,0.1853,0.121,0.1835,0.1655,0.1247,0.0631,0.0477]
downs2 = [0.0406,0.0729,0.1422,0.1329,0.209,0.1671,0.1529,0.0976,0.0537]


vstokesfig = plt.figure()
vstokesax = vstokesfig.add_subplot(111)
vstokesax.set_xlabel('Initial z position (H)', fontsize='16')
vstokesax.set_ylabel('Maximum displacement (H)', fontsize='16')
vstokesax.plot(zposthings, ups, label='Upward displacement')
vstokesax.plot(zposthings, downs, label='Downward displacement')
vstokesax.grid()
plt.legend(loc='upper right')

plt.show()
