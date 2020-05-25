# -*- coding: utf-8 -*-
"""
Trying to reproduce results from: "A fast interpretation of self-potential data using
the depth from extreme points method" by Maurizio Fedi and Mahmoud Ahmed Abbas
Also see Conference Paper Abbas for more testing data
"""
import numpy as np
import matplotlib.pyplot as plt
from fatiando.gravmag import transform
from fatiando import gridder, mesher, utils


# sphere model
K = -600 # mV.m^2
theta = np.deg2rad(90) # deg
z0 = 30 # m
x0 = 0 
m = 1.5


x = np.linspace(-100,100,1000)
nlayers= 25
zz = np.linspace(-40,0,25)

# Equation synthetic SP anomalies 
    
Sxz=[]
for hi in zz:
    num= (x-x0)*np.cos(theta) + (hi-z0)*np.sin(theta)
    den = ((x-x0)**2 + (hi-z0)**2)**m
    Sxz = np.concatenate([Sxz,-K*num/den],axis=0)


shape = (len(x),nlayers)
S_d1x = transform.derivx(x, zz, Sxz, shape, order=1)
S_d1y = transform.derivy(x, zz, Sxz, shape, order=1)

S_d1 = S_d1x - S_d1y
len(S_d1)
len(Sxz)

# ---------------------------------------
# extend the problem to a x,y grid and use  transform.derivz to obtain the first vertical derivative S_d1 ??
# extract then the profile to plot p_S_d1
# ---------------------------------------

# ---------------------------------------
# Try to obtain same result from gradient fct
# 1st order vertical derivative
#dz = np.abs(zz[1] - zz[0])
#S_d1 = np.gradient(Sxz[-1,:],dz)
# ---------------------------------------

# p1 = [-100, 0]
# p2 = [100, 0]
# xx, yy, distance, p_S_d1 = gridder.profile(x, zz, S_d1, p1, p2, 1000)

Sxz = np.reshape(Sxz,[nlayers,len(x)])
S_d1 = np.reshape(S_d1,[nlayers,len(x)])


fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('distance (m)')
ax1.set_ylabel('mV', color=color)
ax1.plot(x, Sxz[-1,:], color=color)
#ax1.plot(x, Sxz, color=color)
ax1.set_ylim([-1,0.2])
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('mV/m', color=color)  # we already handled the x-label with ax1
ax2.plot(x, S_d1[-1,:], color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()

