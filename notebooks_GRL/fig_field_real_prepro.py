# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 18:45:42 2020

@author: Benjamin
simplified from test_field.py
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

# /home/ben/OneDrive/Padova/Software/SourceInversion/Potential_field_imaging/dEXP_imaging/lib/examples/malm/
import notebooks_GRL.load_MALM_model as MALMmod
import notebooks_GRL.load_MALM_real as MALMreal

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


# path =  'E:/Padova/Redaction/Articles/1b_InversionUsingGravityPotMethod/notebooks/data/phNO/'
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
shape = (250,250)
xint_scipy, yint_scipy = gridder.regular((min(Xs)-prl, max(Xs)+prl, 
                          min(Ys)-prl, max(Ys)+prl),shape=shape)




#%% ------------------------------- MALM DATA real
MainPath= r'E:\Padova\Software\SourceInversion\Potential_field_imaging\dEXP_imaging\examples_in_prep\\'
MainPath= '/home/ben/OneDrive/Padova/Software/SourceInversion/Potential_field_imaging/dEXP_imaging/examples_in_prep/'
# os.chdir(MainPath)
## --------- read MALM measured data file (with electrode positions) --------- ##
# RealData = np.loadtxt("./1_Data_2_plot/to_plot.dat",skiprows=0,delimiter='\t') 
out = MALMreal.load_MALM_Porto_real(MainPath + '/malm_models/',
                          MainPath + './malm_models/XYObs_real_f_m3.txt',
                          shape=(100,100),
                          radius=200,
                          rcor=10,
                          rot=0,
                          showfig=True)

# coords_liner = out[3]
p = out[6]         # line points  
p1f , p2f = p
coord_xyz, coord_xyz_int = out[0:2]
xf, yf, zf = coord_xyz
# len(xp)

out_r = MALMreal.load_MALM_Porto_real(MainPath + '/malm_models/',
                          MainPath + './malm_models/XYObs_real_f_m3.txt',
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
Uf = np.copy(Uload_r[0])

# Uf  = np.copy(uf_mirror)
# xf = xf_mirror
# yf = yf_mirror


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


    
XFs, YFs = xint_scipy,yint_scipy

# %%
# prl = 10
# # shape = shape  (max(xp)-min(xp))/
# shape = (30,30)
# xint_scipy, yint_scipy = gridder.regular((min(Xfs)-prl, max(Xfs)+prl, 
#                           min(Yfs)-prl, max(Yfs)+prl),shape=shape)
#%% Solution 1
# extrapolate False and fill with 0 before derivative - mask them later on 
U_int_scipyf = gridder.interp_at(Xfs,Yfs,Uload_r[0], xint_scipy, yint_scipy, algorithm='cubic', extrapolate=False)
InterpData = np.array([xint_scipy, yint_scipy, U_int_scipyf]).T
where_are_NaNs = np.isnan(InterpData)
InterpData[where_are_NaNs] = 0.0074
xint_scipy, yint_scipy, U_int_scipyf_raw = InterpData.T

ax, plt = pEXP.plot_field(XFs,YFs,U_int_scipyf_raw, shape,Vminmax=[0,0.012])
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
ax.set_aspect('equal')
# plt.xlim(min(Xs),max(Ys))
# plt.ylim(min(Xs),max(Ys))
ax.set_xlim(300,500)
ax.set_ylim(300,500)
ax.scatter(Xfs, Yfs, c='k', s=1)
plt.show() 
plt.savefig('publi_raw_field.png', dpi=450)
plt.savefig('publi_raw_field.eps', dpi=450)


#%% Solution 1
# extrapolate False and fill with 0 before derivative - mask them later on 
U_int_scipyf = gridder.interp_at(Xfs,Yfs,Uload_r[1], xint_scipy, yint_scipy, algorithm='cubic', extrapolate=False)
InterpData = np.array([xint_scipy, yint_scipy, U_int_scipyf]).T
where_are_NaNs = np.isnan(InterpData)
InterpData[where_are_NaNs] = 0.0074
xint_scipy, yint_scipy, U_int_scipyf_cor = InterpData.T

ax, plt = pEXP.plot_field(XFs,YFs,U_int_scipyf_cor, shape,Vminmax=[0,0.012])
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
ax.set_aspect('equal')
# plt.xlim(min(Xs),max(Ys))
# plt.ylim(min(Xs),max(Ys))
ax.set_xlim(300,500)
ax.set_ylim(300,500)
ax.scatter(Xfs, Yfs, c='k', s=1)
# plt.show()
# plt.savefig('publi_cor_field.png', dpi=450)
plt.show() 
plt.savefig('publi_cor_field.eps', dpi=450)
plt.savefig('publi_cor_field.svg', dpi=450)

#%% Smoothing 2d

UF = U_int_scipyf_cor
UFs = dEXP.smooth2d(XFs, YFs, UF, sigma=2)
# plt.savefig('smooth2d' + str(file) + '.png', dpi=450)
#UFs = UF

#%% Solution 1
# extrapolate False and fill with 0 before derivative - mask them later on 

ax, plt = pEXP.plot_field(XFs,YFs,UFs, shape,Vminmax=[0,0.012])
ax.plot(coords_liner_s[2:5,0],coords_liner_s[2:5,1],'k')
ax.set_aspect('equal')

# plt.xlim(min(Xs),max(Ys))
# plt.ylim(min(Xs),max(Ys))
ax.set_xlim(300,500)
ax.set_ylim(300,500)
ax.scatter(Xfs, Yfs, c='k', s=1)
plt.show()
plt.savefig('publi_smooth_field.png', dpi=450)
    
#%% Solution 2

# U_int_scipyf = gridder.interp_at(Xfs,Yfs, Uf, xint_scipy, yint_scipy, algorithm='cubic', extrapolate=True)
# len(U_int_scipyf)
# len(Uf)
