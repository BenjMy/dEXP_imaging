import os

from fatiando.vis import mpl #, myv
from fatiando import gridder, mesher, utils
from fatiando.gravmag import prism, imaging, transform
from fatiando.vis.mpl import square

# my own functions
import dEXP as dEXP
from dEXP import _fit
import plot_dEXP as pEXP

# exemples
import exemples.fwd_mag_sphere as magfwd
import exemples.load_grav_model as grav
import exemples.load_MALM_model as MALM

import set_parameters as para

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15

# icsd functions
from icsd3d.importers.read import load_obs, load_geom



#%% ------------------------------- MALM DATA
# path2files="example_2add_later/Landfill_3d/Ano_0_EA/"  # Real A position without anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_BH_EA/"  # Real A position with big hole anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_SH_EA/"  # Real A position  with small hole anomaly
# Main = 'E:/Padova/Software/SourceInversion/icsd_dev/example_2add_later/Landfill_3d/Ano_1_BH_EA/'
# file = 'OAno_synt'

Main = 'E:/Padova/Experiments/GC_2019_Landfill_Porto_MALM/Test/'
file = 'Ano'
interp = True

dataset = MALM.load_MALM_LandfillPorto(path=Main, 
                                       filename=file,
                                       shape = (300,300),
                                       field=True,
                                       interp = interp)
coord_xyz, coord_xyz_int = dataset[0:2]
U = dataset[2]
coords_liner = dataset[3]
shape, max_elevation = dataset[4:6]
p = dataset[6]         # line points                                       

# set imaging pseudo-inversion parameters                                                                        
parameters = para.set_par(shape=shape,max_elevation=max_elevation)

scaled = parameters[0]
SI = parameters[1]
zp, qorder, nlay = parameters[2:5]
minAlt_ridge, maxAlt_ridge = parameters[5:7]

# xp, yp, zp = coord_xyz_int
xp, yp, zp = coord_xyz
# len(xp)
U = U[0] # U_raw, Ucor, U_int, Ucor_int
p1 , p2 = p

#%% ------------------------------- MAG DATA
# -------------------------------  Model
# xp, yp, zp, U, shape, p1, p2, coord= magfwd.load_mag_synthetic()
# max_elevation=2*max(coord[:,2])
# scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)

#%% ------------------------------- GRAVITY DATA
# -------------------------------  Model
# load_grav_synthetic()
# ga, gza = grav.load_grav_pygimli_cylinder()
# shape = (30,30)
# scaled, SI, zp, qorder, nlay, minAlt_ridge, maxAlt_ridge = para.set_par(shape=shape,max_elevation=max_elevation)

#%% ------------------------------- Plot the data 
pEXP.plot_line(xp, yp, U,p1,p2, interp=interp)

#%% ------------------------------- Pad the edges of grids

xp,yp,U, shape = dEXP.pad_edges(xp,yp,U,shape,pad_type=2)
pEXP.plot_line(xp, yp, U,p1,p2, interp=interp)

#%% ------------------------------- Plot the derivatives
xderiv = transform.derivx(xp, yp, U, shape,order=1)
yderiv = transform.derivy(xp, yp, U, shape,order=1)
zderiv = transform.derivz(xp, yp, U, shape,order=1)

# interp = True
pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=interp)
pEXP.plot_line(xp, yp, xderiv ,p1,p2,title='xderiv',savefig=False, interp=interp)
pEXP.plot_line(xp, yp, zderiv ,p1,p2,title='zderiv',savefig=False, interp=interp)

#%% ------- upward continuation of the field data

mesh, label_prop = dEXP.upwc(xp, yp, zp, U, shape, 
                 zmin=0, zmax=max_elevation, nlayers=nlay, 
                 qorder=qorder)

plt, cmap = pEXP.plot_xy(mesh, label=label_prop)
plt.colorbar(cmap)
        
#%% ------------------------------- ridges identification

# import timeit

dEXP.ridges_minmax_plot(xp, yp, mesh, p1, p2,
                                      label=label_prop,
                                      fix_peak_nb=None,
                                      interp=False) 

dfI,dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,interp=interp,
                                      label=label_prop,fix_peak_nb=None) 

# dfI, dfII, dfIII = dEXP.ridges_minmax(xp, yp, mesh, p1, p2,
#                                       label=label_prop,
#                                       # minAlt_ridge=minAlt_ridge,
#                                       maxAlt_ridge=maxAlt_ridge,
#                                       fix_peak_nb=None) 
 
#%% ------------------------------- plot ridges over continuated section

fig = plt.figure()
ax = plt.gca()
pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI,dfII,dfIII,ax=ax)

#%% ------------------------------- filter ridges regionally constrainsted)
   

# dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
#                                            minAlt_ridge,maxAlt_ridge,
#                                            minlength=5,rmvNaN=True)

dfI_f,dfII_f, dfIII_f = dEXP.filter_ridges(dfI,dfII,dfIII,
                                           minAlt_ridge,maxAlt_ridge,
                                           minlength=5,rmvNaN=True)

df_f = dfI_f, dfII_f, dfIII_f

#%% ------------------------------- plot ridges fitted over continuated section

fig = plt.figure()
ax = plt.gca()

pEXP.plot_xy(mesh, label=label_prop, ax=ax) #, ldg=)
pEXP.plot_ridges_harmonic(dfI_f,dfII_f,dfIII_f,ax=ax,label=True)

df_fit = dEXP.fit_ridges(df_f) # fit ridges on filtered data

pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-max_elevation*2,
                          ridge_type=[0,1,2],ridge_nb=None)

# pEXP.plot_ridges_sources(df_fit, ax=ax, z_max_source=-30, 
#                           ridge_type=[0]) # ridge_nb = [[1,2],[1,3]]

#%% ------------------------------- save intersection
# TO IMPLEMENT !
dEXP.ridges_intersection_Z0(df_fit,ax=ax,ridge_nb=[3,4]) # 

# f= fit
# idx = np.argwhere(np.diff(np.sign(f - g))).flatten()
# plt.plot(x[idx], f[idx], 'ro')

# import itertools  
    
# l = ['Geeks', 'for', 'Geeks']  

# import itertools


# # defining iterator  
# iterators = itertools(l)  
    
# # for in loop  
# for i in range(6):  
        
#     # Using next function  
#     print(next(iterators), end = " ")  

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
    