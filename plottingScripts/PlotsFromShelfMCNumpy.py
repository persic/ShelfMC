import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
import matplotlib.patches as mpatches

infn = sys.argv[1]

Array = np.load(infn)

X = Array['xyz'][:,0]
Y = Array['xyz'][:,1]
Z = Array['xyz'][:,2]
polZ = Array['n_pol'][:,2]
Rho = np.sqrt(X*X+Y*Y)
weight = Array['weight']
Exp = Array['log10E']
nEntries = len(Array)

phiHat = np.column_stack((-np.sin(np.deg2rad(Array['arrivalPhi'])),np.cos(np.deg2rad(Array['arrivalPhi'])),np.zeros(nEntries)))
thetaHat = np.column_stack((np.cos(np.deg2rad(Array['arrivalPhi']))*np.cos(np.deg2rad(Array['arrivalTheta'])),np.sin(np.deg2rad(Array['arrivalPhi']))*np.cos(np.deg2rad(Array['arrivalTheta'])),-np.sin(np.deg2rad(Array['arrivalTheta']))))

pol_phi = (Array['n_pol']*phiHat).sum(1)
pol_theta = (Array['n_pol']*thetaHat).sum(1)

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

fig7, ax7 = plt.subplots()
fig7.suptitle('Polz vs Emission Angle')
cb = ax7.hist2d(polZ,Array['viewangle']-Array['coneAngle'],weights=weight,bins=60)
ax7.set_xlabel('cos(theta_pol)')
ax7.set_ylabel('viewaAngle-coneAngle')
plt.colorbar(cb[3],ax=ax7)

fig8, ax8 = plt.subplots()
fig8.suptitle('Pol_theta vs Emission Angle')
cb = ax8.hist2d(np.abs(pol_theta),Array['viewangle']-Array['coneAngle'],weights=weight,bins=60)
ax8.set_xlabel('pol_thetaHat')
ax8.set_ylabel('viewaAngle-coneAngle')
plt.colorbar(cb[3],ax=ax8)

fig9, ax9 = plt.subplots(2,2)
fig9.suptitle('Polz vs Emission Angle')

minTheta = 60.0
maxTheta = 65.0
ThetaFilter = np.all([Array['arrivalTheta']>minTheta,Array['arrivalTheta']<maxTheta],axis=0)
ax9[0][0].set_title('{} - {}deg'.format(minTheta,maxTheta))
cb = ax9[0][0].hist2d(polZ[ThetaFilter],Array['viewangle'][ThetaFilter]-Array['coneAngle'][ThetaFilter],weights=weight[ThetaFilter],bins=60)
ax9[0][0].set_xlabel('cos(theta_pol)')
ax9[0][0].set_ylabel('viewaAngle-coneAngle')
plt.colorbar(cb[3],ax=ax9[0][0])

minTheta = 40.0
maxTheta = 45.0
ThetaFilter = np.all([Array['arrivalTheta']>minTheta,Array['arrivalTheta']<maxTheta],axis=0)
ax9[0][1].set_title('{} - {}deg'.format(minTheta,maxTheta))
cb = ax9[0][1].hist2d(polZ[ThetaFilter],Array['viewangle'][ThetaFilter]-Array['coneAngle'][ThetaFilter],weights=weight[ThetaFilter],bins=60)
ax9[0][1].set_xlabel('cos(theta_pol)')
ax9[0][1].set_ylabel('viewaAngle-coneAngle')
plt.colorbar(cb[3],ax=ax9[0][1])

minTheta = 20.0
maxTheta = 25.0
ThetaFilter = np.all([Array['arrivalTheta']>minTheta,Array['arrivalTheta']<maxTheta],axis=0)
ax9[1][0].set_title('{} - {}deg'.format(minTheta,maxTheta))
cb = ax9[1][0].hist2d(polZ[ThetaFilter],Array['viewangle'][ThetaFilter]-Array['coneAngle'][ThetaFilter],weights=weight[ThetaFilter],bins=60)
ax9[1][0].set_xlabel('cos(theta_pol)')
ax9[1][0].set_ylabel('viewaAngle-coneAngle')
plt.colorbar(cb[3],ax=ax9[1][0])

minTheta = 70.0
maxTheta = 75.0
ThetaFilter = np.all([Array['arrivalTheta']>minTheta,Array['arrivalTheta']<maxTheta],axis=0)
ax9[1][1].set_title('{} - {}deg'.format(minTheta,maxTheta))
cb = ax9[1][1].hist2d(polZ[ThetaFilter],Array['viewangle'][ThetaFilter]-Array['coneAngle'][ThetaFilter],weights=weight[ThetaFilter],bins=60)
ax9[1][1].set_xlabel('cos(theta_pol)')
ax9[1][1].set_ylabel('viewaAngle-coneAngle')
plt.colorbar(cb[3],ax=ax9[1][1])

thetaRanges = [(20,21),(40,41),(60,61)]
nRanges = len(thetaRanges)
rgba = np.zeros((nRanges,nEntries,4))
for i in xrange(nRanges):
    rgba[i,:,0:3]=c.to_rgb('C{}'.format(i))
    rgba[i,:,3]=np.clip(weight,a_min=None,a_max=1.0)

fig10, ax10 = plt.subplots()
ax9[1][1].set_xlabel('cos(theta_pol)')
ax9[1][1].set_ylabel('viewaAngle-coneAngle')

for i,(minTheta,maxTheta) in enumerate(thetaRanges):
    ThetaFilter = np.all([Array['arrivalTheta']>minTheta,Array['arrivalTheta']<maxTheta],axis=0)
    ax10.scatter(polZ[ThetaFilter],Array['viewangle'][ThetaFilter]-Array['coneAngle'][ThetaFilter],color=rgba[i][ThetaFilter])

ax10handles = [mpatches.Patch(color='C{}'.format(i), label='{} - {} deg'.format(thetaRanges[i][0],thetaRanges[i][1])) for i in xrange(nRanges)]
ax10.legend(handles=ax10handles)

#fig7, ax7 = plt.subplots()
#fig7.suptitle('Polz vs Emission Angle')
#ax7.scatter(polZ,Array['viewangle']-Array['coneAngle'],alpha=None,c=Array['arrivalTheta'])
#ax7.set_xlabel('cos(theta_pol)')
#ax7.set_ylabel('viewaAngle-coneAngle')

plt.show()
