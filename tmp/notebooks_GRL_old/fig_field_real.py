# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 09:58:24 2021

@author: Benjamin
"""
from fatiando.gravmag import transform
import pickle

# my own functions
import lib.dEXP as dEXP
import lib.plot_dEXP as pEXP
import lib.utils_dEXP as uEXP
import matplotlib.pyplot as plt
from mpl_axes_aligner import align

plt.rcParams['font.size'] = 15
import numpy as np

#%% ------------------------------- MALM DATA
data_dir= './data/'
#data_dir= 'E:/Padova/Software/SourceInversion/Potential_field_imaging/dEXP_imaging/examples_in_prep/'

# data are saved from test_synth.py
savename = "field_real"
file = open(data_dir+'fig3_data.pkl','rb')
u = pickle._Unpickler(file)
u.encoding = 'latin1'
data = u.load()

XFs, YFs, UF, UFs = data['XYU']
p1_s, p2_s = data['p12']
shape = data['shape']
coords_liner_s = data['coords_liner']
smooth = 'CubicSmoothingSpline' #'Hanning+Lowpass'
# smooth = 'CubicSmoothingSpline + interp1d' #'Hanning+Lowpass'

savename = 'field_real'
interp = False
interp_size = 300
x_axis = 'x'
zp = 0
max_elevation = 30
minAlt_ridge, maxAlt_ridge = 12.5,37.5
nlay = 25
qorder = 0

UF = UFs

#%% ------------------------------- Plot the derivatives

xderiv = transform.derivx(XFs, YFs, UF, shape,order=0)
yderiv = transform.derivy(XFs, YFs, UF, shape,order=0)
zderiv = transform.derivz(XFs, YFs, UF, shape,order=0)

# interp = True
pEXP.plot_line(XFs, YFs, xderiv ,p1_s,p2_s,title='xderiv',x_resolution= interp_size,
            savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

plt.savefig('xderiv' + savename + '.png', dpi=450)

pEXP.plot_line(XFs, YFs, yderiv ,p1_s,p2_s,title='yderiv',x_resolution= interp_size,
            savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

plt.savefig('yderiv' + savename + '.png', dpi=450)

pEXP.plot_line(XFs, YFs, zderiv ,p1_s,p2_s,title='zderiv',x_resolution= interp_size,
            savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

plt.savefig('zderiv' + savename + '.png', dpi=450)

#%% ------- upward continuation of the field data
p = [p1_s,p2_s]

mesh, label_prop = dEXP.upwc(XFs, YFs, zp, UF, shape, 
              zmin=0, zmax=max_elevation, nlayers=nlay, 
              qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis,
                     Vminmax=[0,0.0125], p1p2=p)
cbar = plt.colorbar(cmap, shrink= 0.5, orientation='vertical')
cbar.set_label('upwc voltage ($V.m^2$)')
plt.tight_layout()
plt.savefig('upwc voltage' + savename + '.png', dpi=450)

#%% DEXP ratio
# x_axis = 'x'
qratio = [1,0]
mesh_dexp, label_dexp = dEXP.dEXP_ratio(XFs, YFs, zp, UF, shape, 
              zmin=0, zmax=max_elevation, nlayers=nlay, 
              qorders=qratio)
fig = plt.figure(figsize=(15,3))
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
          markerMax=True,qratio=str(qratio),Vminmax=[0,0.075],
          p1p2=np.array([p1_s,p2_s]), ax=ax, Xaxis=x_axis) #, ldg=)
# plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
#              markerMax=True,qratio=str(qratio)
#              ax=ax, Xaxis=x_axis) #, ldg=)
cbar = plt.colorbar(cmap, shrink= 0.5, orientation='vertical')
cbar.set_label('ratio voltage (V)')

# if x_axis=='y':
# square([xA_r_new[0], xA_r_new[1], -z1, -z2])
# else:   
# square([yA_r[0], yA_r[1], -z1, -z2])
plt.xlim([200,600])
plt.savefig('ratiosfield.png', dpi=450)
# plt.legend('')


# %% ridges identification

dEXP.ridges_minmax_plot(XFs, YFs, mesh, p1_s, p2_s,
                                  label=label_prop,
                                  interp=interp,x_resolution= interp_size,
                                  smooth=smooth,
                                  fix_peak_nb=2,
                                  method_peak='find_peaks',
                                  showfig=True,
                                  Xaxis=x_axis)  
#%%
# or  find_peaks or peakdet or spline_roots
# dfI,dfII, dfIII = dEXP.ridges_minmax(Xs, Ys, mesh, p1_s, p2_s,interp=interp,x_resolution= interp_size,
#                                       label=label_prop,fix_peak_nb=2,
#                                       smooth=smooth, # true = low pass, otherwise specify the filter to apply
#                                       method_peak='find_peaks',
#                                       showfig=True,
#                                       Xaxis=x_axis) 

D = dEXP.ridges_minmax(XFs, YFs, mesh, p1_s, p2_s,
                                  label=label_prop,
                                  method_peak='find_peaks',
                                  fix_peak_nb=2,
                                  returnAmp=True,
                                  showfig=True,
                                  Xaxis=x_axis,
                                  interp=interp,x_resolution= interp_size,
                                  smooth = smooth,
                                  qorder=qorder)  

dfI, dfII, dfIII =  D[0:3]
hI, hII, hIII  = D[3:6]
heights  = D[3:6]

#%% ------------------------------- plot ridges over continuated section

fig = plt.figure(figsize=(15,3))
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis,
        Vminmax=[0,0.0125], p1p2=p)
cbar = plt.colorbar(cmap,shrink= 0.5)
cbar.set_label('upwc voltage (V)')
plt.tight_layout()
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
plt.xlim([200,600])
plt.savefig('ridges_raw_field.png', dpi=450)

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
          Vminmax=[0,0.0125], p1p2=p)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax1,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

ax2 = pEXP.plot_ridges_sources(df_fit, ax=ax1, z_max_source=-max_elevation,
                      ridge_type=[0,1,2],ridge_nb=None)

labels_ax = ax.get_yticks() 
labels_ax= labels_ax[labels_ax>0]

labels_ax2 = ax2.get_yticks() 
labels_ax2= labels_ax2[labels_ax2<0]

ax.set_yticks(labels_ax)
ax2.set_yticks(labels_ax2)

# Adjust the plotting range of two y axes
org1 = 0.0  # Origin of first axis
org2 = 0.0  # Origin of second axis
pos = 0.2  # Position the two origins are aligned
align.yaxes(ax1, org1, ax2, org2, pos)


cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('upwc voltage (V)')
plt.tight_layout()
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax1)
plt.xlim([200,600])


plt.savefig('ridgesfield.png', dpi=450)
