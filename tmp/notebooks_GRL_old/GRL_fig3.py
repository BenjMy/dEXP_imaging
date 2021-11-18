# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:45:42 2020

@author: Benjamin
"""


import os
import pickle

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import lib.dEXP as dEXP
from lib.dEXP import _fit
import lib.plot_dEXP as pEXP
import lib.utils_dEXP as uEXP

import notebooks_GRL.load_MALM_model as MALMmod
import notebooks_GRL.load_MALM_real as MALMreal
from mpl_axes_aligner import align

import lib.set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

# icsd functions
#from importers.read import load_obs, load_geom


#%% ------------------------------- MALM DATA synthetic

interp_size = 300
smooth = 'CubicSmoothingSpline' #'Hanning+Lowpass'
# smooth = 'CubicSmoothingSpline + interp1d' #'Hanning+Lowpass'
interp = False
file = 'NoAno'


path =  './data/phNO/'
dataset = MALMmod.load_MALM_Landfill_model(path=path, 
                                    filename=file,
                                    shape = (300,300),
                                    field=True,
                                    interp = interp,
                                    radius=60) # length of p1p2 profile

coord_xyz, coord_xyz_int = dataset[0:2]
coord_xyz_int
Uload = dataset[2]
coords_liner = dataset[3]
shape, max_elevation = dataset[4:6]

dict_data = dataset[7]
dict_data['AnoBool']

xA = (dict_data['HZ'][0][0]+dict_data['HZ'][0][1])/2
x1 = dict_data['HZ'][0][2]
y1 = dict_data['HZ'][0][0]
y2 = dict_data['HZ'][0][1]

z1 = dict_data['HZ'][1]
z2 = z1 - dict_data['HZ'][2]

p = dataset[6]         # line points                                       
# set imaging pseudo-inversion parameters                                                                        
parameters = para.set_par(shape=shape,max_elevation=max_elevation)

scaled = parameters[0]
SI = parameters[1]
zp, qorder, nlay = parameters[2:5]
minAlt_ridge, maxAlt_ridge = parameters[5:7]

max_elevation = 30
# nlay = 50

# xp, yp, zp = coord_xyz_int
xp, yp, zp = coord_xyz
# len(xp)
Uini = Uload[0] # U_raw, Ucor, U_int, Ucor_int
p1 , p2 = p

rot = 60
origin=(max(xp), min(yp))
point_torotate = np.array([xp, yp])
xp_r, yp_r = MALMreal.rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
Xs = xp_r-min(xp_r)
Ys = yp_r-min(yp_r) 
    
prl = 60
# shape = shape  (max(xp)-min(xp))/
shape = (150,150)
xint_scipy, yint_scipy = gridder.regular((min(Xs)-prl, max(Xs)+prl, 
                          min(Ys)-prl, max(Ys)+prl),shape=shape)
    
#%% ------------------------------- MALM DATA real
MainPath= './data/'
# os.chdir(MainPath)
## --------- read MALM measured data file (with electrode positions) --------- ##
# RealData = np.loadtxt("./1_Data_2_plot/to_plot.dat",skiprows=0,delimiter='\t') 
out = MALMreal.load_MALM_Porto_real(MainPath,
                          MainPath + 'XYObs_real_f_m3.txt',
                          shape=(100,100),
                          radius=200,
                          rcor=10,
                          rot=0,
                          showfig=False)

# coords_liner = out[3]
p = out[6]         # line points  
p1f , p2f = p
coord_xyz, coord_xyz_int = out[0:2]
xf, yf, zf = coord_xyz
# len(xp)

out_r = MALMreal.load_MALM_Porto_real(MainPath,
                          MainPath + 'XYObs_real_f_m3.txt',
                          shape=(100,100),
                          radius=200,
                          rcor=10,
                          rot=60,
                          showfig=False)

# coords_liner = out[3]
p = out_r[6]         # line points  
p1f_r , p2f_r = p
coord_xyz_r, coord_xyz_int_r = out_r[0:2]
Uload_r = out_r[2]
uf = Uload_r[1] # U_raw, Ucor, U_int, Ucor_int
# xf_rtmp, yf_rtmp, zf_rtmp = coord_xyz

# plt.figure()
# plt.scatter(xf_rtmp, yf_rtmp, c=uf, cmap='viridis',vmax=0.05)
# plt.colorbar()
# plt.axis('square')
# plt.show()
  

#%% Mirror field against p1p2

# find point position with respect to line equation defined by p1 and p2
U_af, p_af, bool_abovef, U_bf, p_bf = MALMmod.isabove(xf, yf, uf, 
                                  np.array(p1),np.array(p2))
# Compute p1 and p2 line equation ax + by + c = 0
af, bf, cf = MALMmod.slope(p1,p2)

# Mirror points with respect to p1p2 line
Umirrorf, xy_mirrorf = MALMmod.mirrorU_alongLine(U_af,p_af,bool_abovef,af,bf,cf)

U_a_intf = gridder.interp_at(xy_mirrorf[:,0], xy_mirrorf[:,1], Umirrorf, xf, yf, algorithm='cubic', 
                        extrapolate=True)   
# U_mirror_int = np.copy(U_a_int)
U_mirror_intf = np.copy(uf)
U_mirror_intf[np.where(bool_abovef == True)[0]]= U_a_intf[np.where(bool_abovef == True)[0]]

xf_mirror = np.hstack([xy_mirrorf[:,0],p_af[:,0]])
yf_mirror = np.hstack([xy_mirrorf[:,1],p_af[:,1]])
uf_mirror = np.hstack([Umirrorf,U_af])




#%% choose raw or mirrored field

Uf = np.copy(uf)

Uf  = np.copy(uf_mirror)
xf = xf_mirror
yf = yf_mirror


# %% rotate and rescale all

rot = 60
origin=(max(xp), min(yp))
# origin=(max(Xf), min(Yf))
point_torotate = np.array([xp, yp])
xp_r, yp_r = MALMreal.rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)

point_torotate = np.array([xf, yf])
xf_r, yf_r = MALMreal.rotate_60(origin, point_torotate, angle_deg=rot, showfig=False)
Xfs = xf_r-min(xp_r)
Yfs = yf_r-min(yp_r)

plt.figure()
plt.scatter(xf_r, yf_r, c=Uf, cmap='viridis',vmax=0.0175)
plt.colorbar()
plt.axis('square')
plt.show()

# print(Uf[0:2])

point_torotate = np.array([[dict_data['HZ'][0][0],dict_data['HZ'][0][1]],
                          [dict_data['HZ'][0][2],dict_data['HZ'][0][2]]])
xA_r, yA_r =  MALMreal.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False) #58.511
xA_r = xA_r-min(xp_r)
yA_r = yA_r-min(yp_r)
    
point_torotate = np.array(coords_liner).T
coords_linerx, coords_linery =  MALMreal.rotate_60(origin,point_torotate,angle_deg=rot, showfig=False)
coords_linerx = coords_linerx-min(xp_r)
coords_linery = coords_linery-min(yp_r)
coords_liner_s = np.array([coords_linerx, coords_linery]).T
    

plt.figure()
plt.scatter(Xfs, Yfs, c=Uf, cmap='viridis',vmax=0.05)
plt.colorbar()
plt.axis('square')
plt.show()


    

# %%
# prl = 10
# # shape = shape  (max(xp)-min(xp))/
# shape = (30,30)
# xint_scipy, yint_scipy = gridder.regular((min(Xfs)-prl, max(Xfs)+prl, 
#                           min(Yfs)-prl, max(Yfs)+prl),shape=shape)
#%% Solution 1
# extrapolate False and fill with 0 before derivative - mask them later on 
U_int_scipyf = gridder.interp_at(Xfs,Yfs,Uf, xint_scipy, yint_scipy, algorithm='cubic', extrapolate=False)
InterpData = np.array([xint_scipy, yint_scipy, U_int_scipyf]).T
where_are_NaNs = np.isnan(InterpData)
InterpData[where_are_NaNs] = 0.0074
xint_scipy, yint_scipy, U_int_scipyf = InterpData.T



#%% Solution 2

# U_int_scipyf = gridder.interp_at(Xfs,Yfs, Uf, xint_scipy, yint_scipy, algorithm='cubic', extrapolate=True)
# len(U_int_scipyf)
# len(Uf)

#%% Smoothing 2d

XFs, YFs = xint_scipy,yint_scipy
UF = U_int_scipyf
UF = dEXP.smooth2d(XFs, YFs, UF, sigma=2)
# plt.savefig('smooth2d' + str(file) + '.png', dpi=450)
UFs = UF

#%% 

plt.figure()
plt.scatter(XFs, YFs, c=UF, cmap='viridis')
plt.colorbar()
plt.axis('square')
plt.show()


offset = 255 #260
p1_s = np.array([p1f_r[0] -min(xp_r)+offset,p2f_r[0] -min(xp_r)+offset])
p2_s = np.array([p1f_r[1] -min(yp_r),p2f_r[1] -min(yp_r) ])
ax, plt = pEXP.plot_field(XFs,YFs,UF, shape,
            Vminmax=[0,0.02])
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
ax.plot(p1_s,p2_s,'r')

x_axis = 'x'
# xx, yy, distance, profile, ax, plt = pEXP.plot_line(Xs, Ys, U ,p1_s,p2_s, 
#                                                     interp=True, smooth=True, 
#                                                     xaxis = x_axis)
p1_s = np.array([p1f_r[0] -min(xp_r)+offset,p1f_r[1] -min(yp_r)])
p2_s = np.array([p2f_r[0] -min(xp_r)+offset,p2f_r[1] -min(yp_r) +
             abs(p1f_r[1] -min(yp_r) - (yA_r[0]+yA_r[1])/2)
             -abs(p2f_r[1] -min(yp_r) - (yA_r[0]+yA_r[1])/2)])
# p1_s = np.array([p1[0] -min(xp_r)+260,0])
# p2_s = np.array([p2[0] -min(xp_r)+260,800])
xx, yy, distance, profile, ax,plt = pEXP.plot_line(XFs, YFs, UF ,p1_s,p2_s, 
                                        interp=False,
                                        x_resolution = interp_size,
                                        smooth=smooth, 
                                        xaxis = x_axis,
                                        Vminmax=[0,0.2],
                                        limx=[100,650],
                                        limy=[100,650],
                                        showfig=True)
plt.savefig('profile_field.png', dpi=450)

# xA_r_new = [p1_s[0]+xA_r[0]-xA_r[1], p1_s[0]-xA_r[0]+xA_r[1]]  

#%% ------------------------------- plot publi mirror

ax, plt = pEXP.plot_field(XFs,YFs,UF, shape,Vminmax=[0,0.009])
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
plt.axis('square')
# plt.xlim(min(Xs),max(Ys))
# plt.ylim(min(Xs),max(Ys))
plt.xlim(300,500)
plt.ylim(300,500)
plt.savefig('publi_mirror_field.png', dpi=450)
    

#%% SAVE DATA AND PARAMETERS
    
dict_real_data = { "XYU" : [XFs,YFs,UF, UFs], 
            "prl" : 60,
            "shape": shape,
            "coords_liner": coords_liner_s,
            "p12": [p1_s,p2_s]}


afile = open('fig6_data_m' + '.pkl', 'wb')
pickle.dump(dict_real_data, afile)
afile.close()

interp_size
interp
smooth
x_axis
nlay

#%% ------------------------------- Pad the edges of grids

# xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=0) # reflexion=5
# pEXP.plot_line(xp, yp,U,p1,p2, interp=interp)

#%% ------------------------------- Plot the derivatives

xderiv = transform.derivx(XFs, YFs, UF, shape,order=0)
yderiv = transform.derivy(XFs, YFs, UF, shape,order=0)
zderiv = transform.derivz(XFs, YFs, UF, shape,order=0)

# interp = True
pEXP.plot_line(XFs, YFs, xderiv ,p1_s,p2_s,title='xderiv',x_resolution= interp_size,
            savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

plt.savefig('xderiv' + str(file) + '.png', dpi=450)

pEXP.plot_line(XFs, YFs, yderiv ,p1_s,p2_s,title='yderiv',x_resolution= interp_size,
            savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

plt.savefig('yderiv' + str(file) + '.png', dpi=450)

pEXP.plot_line(XFs, YFs, zderiv ,p1_s,p2_s,title='zderiv',x_resolution= interp_size,
            savefig=False, interp=interp, smooth=smooth,  Xaxis=x_axis)

plt.savefig('zderiv' + str(file) + '.png', dpi=450)

#%% ------- upward continuation of the field data
p = [p1_s,p2_s]

mesh, label_prop = dEXP.upwc(XFs, YFs, zp, UF, shape, 
              zmin=0, zmax=max_elevation, nlayers=nlay, 
              qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop,Xaxis=x_axis,
                     Vminmax=[0,0.0125], p1p2=p)
cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('upwc voltage (V)')
plt.tight_layout()
plt.savefig('upwc voltage' + str(file) + '.png', dpi=450)

#%% DEXP ratio

# x_axis = 'x'
qratio = [1,0]
mesh_dexp, label_dexp = dEXP.dEXP_ratio(XFs, YFs, zp, UF, shape, 
              zmin=0, zmax=max_elevation, nlayers=nlay, 
              qorders=qratio)
fig = plt.figure(figsize=(6, 2), dpi=450)
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
          markerMax=True,qratio=str(qratio),Vminmax=[0,0.075],
          p1p2=np.array([p1_s,p2_s]), ax=ax, Xaxis=x_axis) #, ldg=)
# plt, cmap = pEXP.plot_xy(mesh_dexp, label=label_dexp,
#              markerMax=True,qratio=str(qratio)
#              ax=ax, Xaxis=x_axis) #, ldg=)
cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('ratio\nvoltage (V)')

# if x_axis=='y':
# square([xA_r_new[0], xA_r_new[1], -z1, -z2])
# else:   
# square([yA_r[0], yA_r[1], -z1, -z2])
plt.xlim([200,600])
plt.savefig('fig3b.png', dpi=450)
plt.savefig('fig3b.pdf', dpi=450)
plt.savefig('fig3b.svg', dpi=450)
plt.savefig('fig3b.eps', dpi=450)
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

fig = plt.figure(figsize=(6, 2), dpi=450)
ax = plt.gca()
plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis,
        Vminmax=[0,0.0125], p1p2=p)
cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('upwc\nvoltage (V)')
plt.tight_layout()
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
plt.xlim([200,600])
ax.set_xlabel('y (m)')
ax.set_ylabel('elevation\n(m)')
plt.savefig('fig3a.png', dpi=450)
plt.savefig('fig3a.pdf', dpi=450)
plt.savefig('fig3a.svg', dpi=450)
plt.savefig('fig3a.eps', dpi=450)
# plt.savefig('ridges_raw_field.png', dpi=450)

#%% ------------------------------- filter ridges regionally constrainsted)


dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                        minDepth=minAlt_ridge, 
                                        maxDepth=maxAlt_ridge,
                                        minlength=5,rmvNaN=True,
                                        xmin=100, xmax=700,
                                        Xaxis=x_axis)

df_f = dfI_f, dfII_f, dfIII_f

#%% ------------------------------- plot ridges fitted over continuated section

fig, ax1 = plt.subplots(figsize=(15,3))

plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax, Xaxis=x_axis,
          Vminmax=[0,0.0125], p1p2=p)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=True) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                      ridge_type=[0,1,2],ridge_nb=None)

cbar = plt.colorbar(cmap,shrink=0.25, pad=0.04)
cbar.set_label('upwc voltage (V)')
plt.tight_layout()
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)
plt.xlim([200,600])


plt.savefig('ridgesfield.png', dpi=450)
# square([y1, y2, z1, z2])
#%%

fig, ax1 = plt.subplots(figsize=(15,3))

plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=ax1, Xaxis=x_axis,
          Vminmax=[0,0.0125], p1p2=p)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax1,label=True)

df_fit = dEXP.fit_ridges(df_f, rmvOutliers=False) # fit ridges on filtered data

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