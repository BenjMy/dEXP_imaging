# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:55:57 2020

@author: Benjamin
"""

import matplotlib.pyplot as plt
import numpy as np
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP
import set_parameters as para
import examples.magnetic.fwdmag.fwd_mag_sphere as magfwd


#%%
# Create a model using geometric objects from fatiando.mesher
Atest=[10e3,10e3,2e3]
Atest=[5e3,5e3,2e3]

x_axis = 'y'
rtest=1.5e3
xp, yp, zp, U, shape, p1, p2, coord= magfwd.load_mag_synthetic_ploxy_test(A=Atest, radius=rtest)
max_elevation=2*max(coord[:,2])
scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)
interp = True

#%%
# Plot field data over a 2d line crossing the anomalies
# pEXP.plot_line(xp, yp, U,p1,p2, interp=False, Xaxis='x')

pEXP.plot_line(xp, yp, U,p1,p2, interp=True, Xaxis=x_axis)

#%%
# Upward continuation of the field data with discretisation in altitude controlled by the number of layers (nlay) and the maximum elevation desired (max_elevation)
mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
plt.colorbar(cmap)

p = np.array([p1, p2])
plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=p) #, ldg=)
plt.colorbar(cmap)

# plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis='y', p1p2=p) #, ldg=)
# plt.colorbar(cmap)


# plt, cmap = pEXP.slice_mesh(xp, yp, mesh, label_prop, p1, p2, interp=True, Xaxis='x')
# plt.colorbar(cmap)


#%%
# Ridges identification: plot all extremas obtained via find_peaks function (numpy) for a given altitude
dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks',
                                      showfig=True,
                                      Xaxis=x_axis)  

#%%
# Ridges identification at all levels: plot extremas obtained via find_peaks function (numpy) for all 3 types of extremas familly RI, RII and RIII
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      method_peak='find_peaks',
                                      showfig=True,
                                      Xaxis=x_axis)  

#%%
# filter ridges using a minimum length criterium and and filter for a specific range of altitude
dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            1,maxAlt_ridge,
                                            minlength=8,rmvNaN=True)
df_f = dfI_f, dfII_f, dfIII_f

#%%
# plot filtered ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

plt, cmap = pEXP.plot_xy(mesh, label=label_prop, Xaxis=x_axis, p1p2=p, ax=ax) #, ldg=)
plt.colorbar(cmap)


pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                          ridge_type=[0,1,2],ridge_nb=None)

