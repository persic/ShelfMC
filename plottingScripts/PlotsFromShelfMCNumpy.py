import sys, os
import numpy as np
import matplotlib.pyplot as plt

infn = sys.argv[1]

Array = np.load(infn)

X = Array['xyz'][:,0]
Y = Array['xyz'][:,1]
Z = Array['xyz'][:,2]
polZ = Array['n_pol'][:,2]
Rho = np.sqrt(X*X+Y*Y)
weight = Array['weight']
Exp = Array['log10E']

xBins = np.arange(0,4050,50)
yBins = np.arange(0,2750,50)

fig1, ax1 = plt.subplots()
fig1.suptitle('Vertex Location')
sc = ax1.hist2d(Rho,2700.0-Z,weights=weight,bins=[xBins,yBins])
ax1.set_xlabel('Horizontal Distance (m)')
ax1.set_ylabel('Depth (m)')
ax1.set_ylim([2700,0])
plt.colorbar(sc[3],ax=ax1)

fig2, ax2 = plt.subplots()
fig2.suptitle('Arrival Zenith')
ax2.set_xlabel('Theta (deg)')
ax2.hist(Array['arrivalTheta'],bins=np.arange(0,93,3))

fig3, ax3 = plt.subplots()
fig3.suptitle('Arrival Azimuth')
ax3.set_xlabel('Phi (deg)')
ax3.hist(Array['arrivalPhi'],bins=np.arange(0,370,10))

fig4, ax4 = plt.subplots()
fig4.suptitle('Voltage')
ax4.hist(Array['vMaxTrig']*1e6,weights=weight,bins=60)
ax4.set_xlabel('MaxV (micro-volt)')

fig5, ax5 = plt.subplots()
fig5.suptitle('Emission Angle')
ax5.hist(Array['viewangle']-Array['coneAngle'],weights=weight,bins=60)
ax5.set_xlabel('viewaAngle-coneAngle')

fig6, ax6 = plt.subplots()
fig6.suptitle('Vertical Polarization')
ax6.hist(polZ,weights=weight,bins=np.linspace(0,1,60))
ax6.set_xlabel('cos(theta_pol)')

plt.show()
