# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:45:42 2020

@author: Benjamin
"""

import pickle
from fatiando.vis.mpl import square
from fatiando.gravmag import transform

# my own functions
import lib.dEXP as dEXP
import lib.plot_dEXP as pEXP

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

from mpl_axes_aligner import align

# icsd functions
#from icsd3d.importers.read import load_obs, load_geom


#%% ------------------------------- MALM DATA
data_dir= './data/'
#data_dir= 'E:/Padova/Software/SourceInversion/Potential_field_imaging/dEXP_imaging/examples_in_prep/'

# data are saved from test_synth.py
savename = "field_model"

file2load='fig5_data'
# NoAno_synth_landflill_data
# Ano_synth_landflill_data
# fig5_data
# NoAno_fig5_data
# Ano_fig5_data

file = open(data_dir+file2load+'.pkl','rb')
u = pickle._Unpickler(file)
u.encoding = 'latin1'
data = u.load()

Xs, Ys, U = data['XYU']
p1_s, p2_s = data['p12']
shape = data['shape']
coords_liner_s = data['coords_liner']
smooth = 'CubicSmoothingSpline' #'Hanning+Lowpass'
# smooth = 'CubicSmoothingSpline + interp1d' #'Hanning+Lowpass'
xA_r_new, yA_r, z1, z2 = data['xAyAzA1zA2']


interp = False
interp_size = 300
x_axis = 'x'
zp = 0
max_elevation = 30
minAlt_ridge, maxAlt_ridge = 12.5,37.5
nlay = 25
qorder = 0

#%% plot potential field

xx, yy, distance, profile, ax,plt = pEXP.plot_line(Xs, Ys, U ,p1_s,p2_s, 
                                            interp=False,
                                            x_resolution = interp_size,
                                            smooth=smooth, 
                                            xaxis = x_axis,
                                            Vminmax=[0,0.35],
                                            limx=[100,650],
                                            limy=[100,650],
                                            showfig=True)
plt.savefig('profile' + savename + '.png', dpi=450)


#%% ------------------------------- plot publi mirror

ax, plt = pEXP.plot_field(Xs,Ys,U, shape,Vminmax=[0,0.35])
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
plt.axis('square')
# plt.xlim(min(Xs),max(Ys))
# plt.ylim(min(Xs),max(Ys))
plt.xlim(300,500)
plt.ylim(300,500)
plt.savefig('publi_mirror' + savename + '.png', dpi=450)
plt.savefig('fig2c' + savename + '.png', dpi=450)
plt.savefig('fig2c.svg', dpi=450)
plt.savefig('fig2c.pdf', dpi=450)
plt.savefig('fig2c.png', dpi=450)

#%% ------------------------------- Pad the edges of grids

# xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
# pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)

#%% ------------------------------- Plot the derivatives

xderiv = transform.derivx(Xs, Ys, U, shape,order=0)
yderiv = transform.derivy(Xs, Ys, U, shape,order=0)
zderiv = transform.derivz(Xs, Ys, U, shape,order=0)

# interp = True
#pEXP.plot_line(Xs, Ys, xderiv ,p1_s,p2_s,title='xderiv',x_resolution= interp_size,
#                savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

#plt.savefig('xderiv' + savename + '.png', dpi=450)

#pEXP.plot_line(Xs, Ys, yderiv ,p1_s,p2_s,title='yderiv',x_resolution= interp_size,
#                savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

#plt.savefig('yderiv' + savename + '.png', dpi=450)

#pEXP.plot_line(Xs, Ys, zderiv ,p1_s,p2_s,title='zderiv',x_resolution= interp_size,
#                savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

#plt.savefig('zderiv' + savename + '.png', dpi=450)

#%% ------- upward continuation of the field data
p = [p1_s,p2_s]

mesh, label_prop = dEXP.upwc(Xs, Ys, zp, U, shape, 
                  zmin=0, zmax=max_elevation, nlayers=nlay, 
                  qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis,
                         Vminmax=[0,0.35], p1p2=p)

cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('upwc voltage (V)')
plt.tight_layout()
plt.savefig('upwc voltage' + savename + '.png', dpi=450)

#%% DEXP ratio
# x_axis = 'x'
qratio = [1,0]
mesh_dexp, label_dexp = dEXP.dEXP_ratio(Xs, Ys, zp, U, shape, 
                  zmin=0, zmax=max_elevation, nlayers=nlay, 
                  qorders=qratio)
fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
              markerMax=True,qratio=str(qratio),Vminmax=[0,0.075],
              p1p2=np.array([p1_s,p2_s]), ax=ax, Xaxis=x_axis) #, ldg=)
# plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
#              markerMax=True,qratio=str(qratio)
#              ax=ax, Xaxis=x_axis) #, ldg=)
cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('ratio voltage (V)')

if x_axis=='y':
    square([xA_r_new[0], xA_r_new[1], -z1, -z2])
else:   
    square([yA_r[0], yA_r[1], -z1, -z2])
plt.xlim([200,600])
plt.savefig('ratios_' + savename + '.png', dpi=450)



# %% ridges identification

dEXP.ridges_minmax_plot(Xs, Ys, mesh, p1_s, p2_s,
                                      label=label_prop,
                                      interp=interp,x_resolution= interp_size,
                                      smooth=smooth,
                                      fix_peak_nb=2,
                                      method_peak='find_peaks',
                                      showfig=False,
                                      Xaxis=x_axis)  
#%%

D = dEXP.ridges_minmax(Xs, Ys, mesh, p1_s, p2_s,
                                      label=label_prop,
                                      method_peak='find_peaks',
                                      fix_peak_nb=2,
                                      returnAmp=True,
                                      showfig=False,
                                      Xaxis=x_axis,
                                      interp=interp,x_resolution= interp_size,
                                      smooth = smooth,
                                      qorder=qorder)  

dfI, dfII, dfIII =  D[0:3]
hI, hII, hIII  = D[3:6]
heights  = D[3:6]

#%% ------------------------------- plot ridges over continuated section
    
fig = plt.figure()
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis,
            Vminmax=[0,0.35], p1p2=p)
cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('upwc voltage (V)')
plt.tight_layout()
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
plt.xlim([200,600])
plt.savefig('ridges_raw_' + savename + '.png', dpi=450)

#%% ------------------------------- filter ridges regionally constrainsted)


dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                            minDepth=minAlt_ridge, maxDepth=maxAlt_ridge,
                                            minlength=5,rmvNaN=True,
                                            xmin=100, xmax=700,
                                            Xaxis=x_axis)

df_f = dfI_f, dfII_f, dfIII_f

#%% ------------------------------- plot ridges fitted over continuated section

fig, ax1 = plt.subplots(figsize=(15,3))

plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax1, Xaxis=x_axis,
              Vminmax=[0,0.35], p1p2=p)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax1,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

ax2 = pEXP.plot_ridges_sources(df_fit, ax=ax1, z_max_source=-max_elevation*2,
                          ridge_type=[0,1,2],ridge_nb=None)

cbar = plt.colorbar(cmap,shrink=0.7)
cbar.set_label('upwc voltage ($V.m^2$)')
plt.tight_layout()
#pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax1)
plt.xlim([200,600])
ax1.set_ylim([0,30])
ax2.set_ylim([-30,0])

labels_ax1 = ax1.get_yticks() 
labels_ax1= labels_ax1[labels_ax1>0]

labels_ax2 = ax2.get_yticks() 
labels_ax2= labels_ax2[labels_ax2<0]

ax1.set_yticks(labels_ax1)
ax2.set_yticks(labels_ax2)

# Adjust the plotting range of two y axes
org1 = 0.0  # Origin of first axis
org2 = 0.0  # Origin of second axis
pos = 0.5  # Position the two origins are aligned
align.yaxes(ax1, org1, ax2, org2, pos)



#ax2.spines['right'].set_position(('axes', 400))

if x_axis=='y':
    square([xA_r_new[0], xA_r_new[1], z1, z2])
else:   
    square([yA_r[0], yA_r[1], z1, z2])



plt.savefig('ridges_' + savename + '.png', dpi=450)
