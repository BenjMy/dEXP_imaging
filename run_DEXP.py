import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

import dEXP as dEXP
import plot_dEXP as pEXP
import exemples.fwd_mag_sphere as magfwd
import exemples.load_MALM_model as MALM

import set_parameters as para

import numpy as np
import matplotlib.pyplot as plt


#%% ------------------------------- MALM DATA

# xp, yp, z, U = MALM.load_MALM()
# scaled, shape, SI, zp, qorder, maxdepth, nlay, minAlt_ridge, maxAlt_ridge = para.set_par()
# p1, p2 = [min(xp), 0], [max(xp), 0]

#%% ------------------------------- MAG DATA
# -------------------------------  Model

A= [10e3,10e3,2e3]
B= [25e3,10e3,1e3]
coord= np.array([A,B])
radius = 1.5e3*np.ones(len(coord))
modelmag = magfwd.anomag_model(coord,radii=radius,inc=50, dec=-30)
xp, yp, zp, field2d, shape = magfwd.fwd_model(modelmag, shape = (300, 300),area = [0, 30e3, 0, 30e3])
magfwd.plot_model(xp, yp, field2d, shape)


U = field2d[1]
p1, p2 = [min(xp), 10e3], [max(xp), 10e3]

scaled, shape, SI, zp, qorder, max_elevation, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=2*max(coord[:,2]))

#%% ------------------------------- GRAVITY DATA
# -------------------------------  Model
# name = '1000_zbot3000_data'
# xp, yp, zp, U= np.loadtxt('./grav_models/' + name + '.txt', unpack=True)
# shape = (25,25) # data interpolation
# maxdepth=10000 # max depth for upward continuation
# minAlt_ridge = 1000
# maxAlt_ridge = 5000
# SI = 2 # structural index
# zp=0  # initial depth (conditionned upward or downward)


#%% ------------------------------- Plot the data
# p1, p2 = [min(xp), 0], [max(xp), 0]
pEXP.plot_line(xp, yp, U,p1,p2)
# #plt.title(strname +'ztop' + str(za) +'_zbot'+ str(zb) + '_data', fontsize=20)
# #plt.savefig(pathFig+strname + '_ExempleFig_z' + str(za) + str(zb) + '_data' + '.png')
# #x1, x2, y1, y2, z1, z2 = np.array(model[0].get_bounds())q
# square([y1, y2, x1, x2])

#%% ------------------------------- Pad the edges of grids

# padtypes = ['0', 'mean', 'edge', 'lintaper', 'reflection', 'oddreflection',
#             'oddreflectiontaper']
# fig = plt.figure()
# ax = plt.gca()

# xs = xp.reshape(shape)
# ys = yp.reshape(shape)
# data = U.reshape(shape)

# padtype = padtypes[3]
# padded_data, nps = gridder.pad_array(data, padtype=padtype)
# # Get coordinate vectors
# pad_x, pad_y = gridder.pad_coords([xs, ys], shape, nps)
# padshape = padded_data.shape
# ax.set_title(padtype)
# ax.pcolormesh(pad_y.reshape(padshape), pad_x.reshape(padshape),
#               padded_data, cmap='RdBu_r')
# ax.set_xlim(pad_y.min(), pad_y.max())
# ax.set_ylim(pad_x.min(), pad_x.max())
    
# shape = padded_data.shape
# U = padded_data.reshape(shape[0]*shape[1])
# xp = pad_x
# yp = pad_y

# p1, p2 = [min(xp), 0], [max(xp), 0]


#%% ------------------------------- Plot the derivatives
xderiv = transform.derivx(xp, yp, U, shape,order=1)
yderiv = transform.derivy(xp, yp, U, shape,order=1)
zderiv = transform.derivz(xp, yp, U, shape,order=1)

ax1 = pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=True)

#%% ------- upward continuation of the field data

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
plt.colorbar(cmap)
        
#%% ------------------------------- ridges identification

import timeit
dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
                                      label=label_prop) 
stop = timeit.default_timer()
# dfI.head(5)
# dfII.head(5)
# dfIII.head(5)

#%% ------------------------------- plot ridges over continuated section

fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)


dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                           minAlt_ridge,maxAlt_ridge,
                                           minlength=5,rmvNaN=True)

df_f = np.array([dfI_f, dfII_f, dfIII_f])

#%% ------------------------------- plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)
points, fit = dEXP.fit_ridges(df_f) # fit ridges


pEXP.plot_ridges_sources(points, fit, ax=ax, z_max_source=-8000,ridge_nb=[5,16,17])
# pEXP.plot_ridges_sources(points, fit, ax=ax, z_max_source=-8000)




#%% ------------------------------- save intersection
# TO IMPLEMENT !
dEXP.ridges_intersection_Z0(fit,ax=ax,ridge_nb=[3,4]) # 

f= fit
idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
plt.plot(x[idx], f[idx], 'ro')

import itertools  
    
l = ['Geeks', 'for', 'Geeks']  

import itertools


# defining iterator  
iterators = itertools(l)  
    
# for in loop  
for i in range(6):  
        
    # Using next function  
    print(next(iterators), end = " ")  

# np.diff(df[0]['EX_xpos1'])
# np.diff(df[0]['EX_xpos2']).any()

# print(np.diff(df[i][k[1]]))

#%% ------------------------------- ridges analysis

# points, fit = dEXP.scalFUN(dfI_f,EXTnb=[1],z0=-10000)
z0est = [0,30,60]
P = []
F = []
for zi in z0est:
    points, fit = dEXP.scalFUN(dfII,EXTnb=[1],z0=zi)
    P.append(points)
    F.append(fit)
    
dfII.isnull().values.any()

# points2, fit = dEXP.scalFUN(dfII_f,EXTnb=[1],z0=3)

# fig = plt.figure()
# ax = plt.gca()
# ax1 = plt.subplot(3,1,1)
# pEXP.plot_scalFUN(P, F, ax=ax1, z0=z0est)



# z0est = [0,30,60]
# P = []
# F = []
# for zi in z0est:
#     points, fit = dEXP.scalEuler(dfII,EXTnb=[1],z0=zi)
#     P.append(points)
#     F.append(fit)
    



#%% ------- dEXP continuation of the field data
mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
                          zmin=0, zmax=max_elevation, nlayers=nlay, 
                          qorder=qorder, SI=1)
# pEXP.plot_z(mesh)
pEXP.plot_xy(mesh, label=label_prop,markerMax=True)      



#%% sum up plot

# fig = plt.figure()
# ax = plt.gca()

# for orderi in range(3):
#     i = orderi +1
#     axupwc = plt.subplot(2, 3, i)
#     mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
#                       zmin=0, zmax=max_elevation, nlayers=nlay, 
#                       qorder=orderi) 
#     plt, cmap = pEXP.plot_xy(mesh, label=label_prop, ax=axupwc)
#     plt.colorbar(cmap)
#     axdexp = plt.subplot(2, 3, i+3)
#     mesh, label_prop = dEXP.dEXP(xp, yp, zp, U, shape, 
#                                   zmin=0, zmax=max_elevation, nlayers=nlay, 
#                                   qorder=orderi, SI=1)
#     pEXP.plot_xy(mesh, label=label_prop,markerMax=True, ax=axdexp) 
#     plt.colorbar(cmap)
    