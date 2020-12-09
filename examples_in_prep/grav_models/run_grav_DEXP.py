"""
Example of gravimetric potential field data analysis using pyDEXP
-----------------------------------------------------------------

This code shows a step-by-step processing of potential field imaging aiming at giving an estimate of the depth of the anomaly depth using the dEXP tranformation method.
dEXP method implementation from Fedi et al. 2012. 
Calculations used :mod:`dEXP`, while plotting use the :mod:`plot_dEXP` module.

The gravimetric model data was created using geometric objects from :mod:`fatiando.mesher`. The forward simulation of the data was done using :mod:`fatiando.gravmag` module.

Sources locations:
    Center of mass = [,,] # xyz coordinates
    l.w.h = ,,, # length, width, height (in m)
Sources properties: 
    density contrast= ?

Implements the DEXP method described in Fedi and Pilkington (2012). Application on a anomaly of density  (gravimetry).

.. note::

    This is part of a larger project aiming at inverting current sources density (see more at: https://icsd-dev.readthedocs.io/en/latest/)


**References**

Uieda, L., V. C. Oliveira Jr, and V. C. F. Barbosa (2013), Modeling the Earth with Fatiando a Terra, Proceedings of the 12th Python in Science Conference, pp. 91 - 98.

Uieda, L. (2018). Verde: Processing and gridding spatial data using Green's functions. Journal of Open Source Software, 3(29), 957. doi:10.21105/joss.00957

Fedi, M., and M. Pilkington (2012), Understanding imaging methods for potential
field data, Geophysics, 77(1), G13, doi:10.1190/geo2011-0078.1

----
"""

import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para

# exemples
import examples.gravimetry.loadgrav.grav_models as grav

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15


#%% 
# load model previously generated using Fatiando a terra package

data_struct = grav.load_grav_fatiando(name='loadgrav/za3000_zb3500_l500_ofs0_dens1200')
xp,yp,zp,U = data_struct['xyzg']
shape = data_struct['shape']
model = data_struct['model']
dens  = data_struct['density']
# scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)

x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())

# p1 =[min(yp),0]
# p2 =[max(yp),0]

p1 =[-6000,0]
p2 =[6000,0]

max_elevation=z2*1.2
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True
qorder = 0

max_elevation/nlay
#%% 
# Plot the data 
xx, yy, distance, profile, ax, plt = pEXP.plot_line(xp, yp, U,p1,p2, interp=True)
zderiv = transform.derivz(xp, yp, U, shape,order=1)

xx, yy, distance, dz, ax, plt = pEXP.plot_line(xp, yp, zderiv, p1,p2, interp=True, title='zderiv')


fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('x(m)')
ax1.set_ylabel('Amplitude of the\n potential field', color=color)
ax1.plot(xx, profile, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('dz', color=color)  # we already handled the x-label with ax1
ax2.plot(xx, dz, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim([-0.0001,3e-4])

fig.tight_layout()  # otherwise the right y-label is slightly clipped
# ax1.set_aspect(aspect=1e4)




#%% 
# Pad the edges of grids (if necessary)

xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
p1 =[min(yp),0]
p2 =[max(yp),0]
x_axis='y'

pEXP.plot_line(xp, yp,U,p1,p2, interp=interp,Xaxis=x_axis)

#%% 
# Upward continuation of the field data

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis)
plt.colorbar(cmap)

#%%
# ridges identification
# dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
#                                       label=label_prop,
#                                       fix_peak_nb=2,
#                                       method_peak='find_peaks')  

# or  find_peaks or peakdet or spline_roots
D = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                        label=label_prop,
                        fix_peak_nb=2,
                        method_peak='find_peaks',
                        returnAmp=True,
                        showfig=True,
                        Xaxis=x_axis)  
dfI, dfII, dfIII =  D[0:3]
hI, hII, hIII  = D[3:6]
H  = D[3:6]

#%% 
# Plot ridges over continuated section

fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)


#%% 
# Filter ridges regionally constrainsted)
   

D_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                            minDepth=1000,
                            maxDepth=3000,
                            minlength=3,rmvNaN=True,
                            heights=[hI, hII, hIII])
    
dfI_f, dfII_f, dfIII_f =  D_f[0:3]
hI_f, hII_f, hIII_f = D_f[3:6]
df_f = D_f[0:3]


#%% 
# Plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

# pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*1.2,
#                           ridge_type=[0,1,2],ridge_nb=None)
pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-6000,
                          ridge_type=[0,1,2],ridge_nb=None)
square([x1, x2, -z1, -z2])
plt.annotate(dens,[(x1 + x2)/2, -(z1+z2)/2])

#%% 
#  ridges analysis: scaling function to determine the SI index

z0 = -(z1 + z2)/2 + 100 # choose an estimate of the depth of the anomaly
# z0 = -100 # choose an estimate of the depth of the anomaly
df_height = D_f[3:6]

ncol = 0
for r_type in range(len(df_height)): # loop over ridges type I, II, III
    ncol = ncol + df_f[r_type].shape[1]-1


fig, axs = plt.subplots(1,ncol, figsize=(15, 6), facecolor='w', edgecolor='k')
# fig.subplots_adjust(hspace = .5, wspace=.001)
axs = axs.ravel()
   
# for i in range(ncol):
# df_f[1]['elevation'].iloc[1]
# dzz = np.log(df_f[1]['elevation'].iloc[1])-np.log(df_f[0]['elevation'].iloc[0])
# dzz = np.log(df_f[1]['elevation'].iloc[1]-df_f[0]['elevation'].iloc[0])

nc = 0
SI_est = []
for r_type in range(len(df_height)): # loop over ridges type I, II, III
    for k in enumerate(df_height[r_type].columns[1:]): # loop over ridges of the same familly
        # df_f[r_type].columns[1:]    
        points, fit, SI_est_tmp , EXTnb = dEXP.scalFUN(df_height[r_type],EXTnb=k[1],z0=z0,
                                                       rmvOutliers=True)
        pEXP.plot_scalFUN(points, fit, z0=z0, 
                          label='R'+ str(r_type+1) + k[1], 
                          ax=axs[nc]) # scaling 
        SI_est.append(SI_est_tmp)
        
        nc = nc + 1
    
SI_mean = np.mean(np.abs(SI_est))

#%% 
# #  ridges analysis
# df_f[1]
# Tau = np.gradient(np.log(up_f_Centralridge)) / np.gradient(np.log(z_r))


#%% 
#  ridges analysis
SI_mean = 3
mesh_dexp, label_dexp = dEXP.dEXP(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder,
                 SI=SI_mean)

fig = plt.figure()
ax = plt.gca()

plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
             markerMax=True,SI=SI_mean,
             p1p2=np.array([p1,p2]), ax=ax, Xaxis=x_axis) #, ldg=)
plt.colorbar(cmap)

if x_axis=='y':
    square([x1, x2, z1, z2])
    plt.annotate(dens,[(x1 + x2)/2, -(z1+z2)/2])
else:   
    square([y1, y2, z1, z2])
    plt.annotate(dens,[(y1 + y2)/2, -(z1+z2)/2])


#%% 

uEXP.multipage(dataname + '_test_DEXP_SI.pdf',  dpi=25)
